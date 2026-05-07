#include "NewtonsMethod.h"

#include <symx>
#include <fmt/format.h>
#include <BlockedSparseMatrix/solve_pcg.h>
#include <Eigen/SparseLU>

#ifdef STARK_HAS_MKL
#include <Eigen/PardisoSupport>
#endif

stark::core::NewtonState stark::core::NewtonsMethod::solve(const double& dt, symx::GlobalPotential& global_energy, symx::SecondOrderCompiledGlobal& compiled_global, Callbacks& callbacks, const Settings& settings, Console& console, Logger& logger)
{
	this->global_energy = &global_energy;
	this->compiled_global = &compiled_global;
	this->callbacks = &callbacks;
	this->settings = &settings;
	this->console = &console;
	this->logger = &logger;

	const int ndofs = global_energy.get_total_n_dofs();

	this->u0.resize(ndofs);
	this->u1.resize(ndofs);
	this->residual.resize(ndofs);

	this->step_newton_it = 0;
	this->step_line_search_count = 0;
	this->cg_iterations_in_step = 0;

	double residual_max = std::numeric_limits<double>::max();
	NewtonState newton_state = NewtonState::Running;
	while (newton_state == NewtonState::Running) {
		this->step_newton_it++;

		if (this->step_newton_it == settings.newton.max_newton_iterations) {
			newton_state = NewtonState::TooManyNewtonIterations;
			break;
		}

		if (settings.debug.symx_finite_difference_check) {
			this->compiled_global->test_derivatives_with_FD(1e-4);
			this->compiled_global->test_derivatives_with_FD(1e-6);
			this->compiled_global->test_derivatives_with_FD(1e-8);
		}

		console.print(fmt::format("\n\t\t {:d}. ", this->step_newton_it), ConsoleVerbosity::NewtonIterations);

		EvalResult result = this->_evaluate_E_grad_hess();

		const double E0 = result.E;

		this->residual = this->_compute_residual(result.grad, dt);
		residual_max = this->residual.maxCoeff();
		console.print(fmt::format("r0 = {:.2e} | ", residual_max), ConsoleVerbosity::NewtonIterations);

		if (residual_max < this->settings->newton.epsilon_residual) {
			newton_state = NewtonState::Successful;
			break;
		}

		const bool linear_system_success = this->_solve_linear_system(this->du, result, dt);
		if (!linear_system_success) {
			newton_state = NewtonState::LinearSystemFailure;
			break;
		}

		const double max_da = this->_compute_acceleration_correction(this->du.cwiseAbs().maxCoeff(), dt);
		console.print(fmt::format("da = {:.2e} | ", max_da), ConsoleVerbosity::NewtonIterations);

		const double du_dot_grad = this->du.dot(result.grad);
		if (du_dot_grad > 0.0) {
			newton_state = NewtonState::LineSearchDoesntDescend;
			break;
		}

		double step_valid_configuration = this->_inplace_max_step_in_search_direction(this->du);
		if (step_valid_configuration < 0.01) {
			callbacks.run_on_intermidiate_state_invalid();
			newton_state = NewtonState::InvalidIntermediateConfiguration;
			break;
		}

		EvalResult result_after = this->_evaluate_E_grad();
		this->residual = this->_compute_residual(result_after.grad, dt);
		residual_max = this->residual.maxCoeff();
		console.print(fmt::format(" r1 = {:.2e}", residual_max), ConsoleVerbosity::NewtonIterations);
		console.print(fmt::format(" energy = {:.2e}", this->_evaluate_E().E), ConsoleVerbosity::NewtonIterations);

		const double du_norm = this->du.array().abs().maxCoeff() * dt;
		console.print(fmt::format(" du_norm = {:.2e}", du_norm), ConsoleVerbosity::NewtonIterations);

		if (settings.newton.use_du_norm_threshhold && settings.newton.use_residual_threshhold) {
			if (du_norm < settings.newton.du_norm_threshhold && residual_max < this->settings->newton.residual.tolerance) {
				newton_state = NewtonState::Successful;
				break;
			}
		} else {
			if (settings.newton.use_du_norm_threshhold && du_norm < settings.newton.du_norm_threshhold) {
				newton_state = NewtonState::Successful;
				break;
			}
			if (settings.newton.use_residual_threshhold && residual_max < this->settings->newton.residual.tolerance) {
				newton_state = NewtonState::Successful;
				break;
			}
		}

		double step_that_worked = this->_inplace_backtracking_line_search(this->du, E0, result_after.E, step_valid_configuration, du_dot_grad);
		if (step_that_worked == 0.0) {
			newton_state = NewtonState::TooManyLineSearchIterations;
			break;
		}
	}

	if (!callbacks.run_is_converged_state_valid()) {
		newton_state = NewtonState::InvalidConvergedState;
	}

	console.print("\n\t\t", ConsoleVerbosity::NewtonIterations);
	console.print(fmt::format("#newton: {:d} ", this->step_newton_it), ConsoleVerbosity::TimeSteps);
	console.print(" | #CG/newton: " + std::to_string((int)(this->cg_iterations_in_step / this->step_newton_it)), ConsoleVerbosity::TimeSteps);
	console.print(" | #line_search/newton: " + std::to_string((int)(this->step_line_search_count / this->step_newton_it)), ConsoleVerbosity::TimeSteps);
	if (newton_state == NewtonState::Successful) {
		console.print(" | converged", ConsoleVerbosity::TimeSteps);
	} else {
		console.print(" | not converged", ConsoleVerbosity::TimeSteps);
	}

	switch (newton_state) {
	case NewtonState::TooManyNewtonIterations:
		console.print(fmt::format("\n\t\t -> Max Newton iterations reached ({:d}) with residual_max {:.2e}. ", settings.newton.max_newton_iterations, residual_max), ConsoleVerbosity::TimeSteps);
		break;
	case NewtonState::TooManyLineSearchIterations:
		console.print(fmt::format("\n\t\t -> Max line search iterations reached ({:d}). ", settings.newton.max_line_search_iterations), ConsoleVerbosity::TimeSteps);
		break;
	case NewtonState::LinearSystemFailure:
		console.print("\n\t\t -> Linear system couldn't find a solution. ", ConsoleVerbosity::TimeSteps);
		break;
	case NewtonState::InvalidIntermediateConfiguration:
		console.print("\n\t\t -> Invalid intermediate configuration couldn't be avoided. ", ConsoleVerbosity::TimeSteps);
		break;
	case NewtonState::LineSearchDoesntDescend:
		console.print("\n\t\t -> Line search doesn't descend. ", ConsoleVerbosity::TimeSteps);
		break;
	case NewtonState::InvalidConvergedState:
		console.print("\n\t\t -> Converged state is not valid. ", ConsoleVerbosity::TimeSteps);
		break;
	default:
		break;
	}
	if (newton_state != NewtonState::Successful) {
		this->console->print_error_msg_and_clear(ConsoleVerbosity::TimeSteps);
	}

	logger.add_to_counter("newton_iterations", this->step_newton_it);
	logger.add_to_counter("CG_iterations", this->cg_iterations_in_step);
	logger.add_to_counter("line_search_iterations", this->step_line_search_count);

	return newton_state;
}

void stark::core::NewtonsMethod::_run_before_evaluation()
{
	this->logger->start_timing("before_energy_evaluation");
	this->callbacks->run_before_energy_evaluation();
	this->logger->stop_timing_add("before_energy_evaluation");
}
void stark::core::NewtonsMethod::_run_after_evaluation()
{
	this->logger->start_timing("after_energy_evaluation");
	this->callbacks->run_after_energy_evaluation();
	this->logger->stop_timing_add("after_energy_evaluation");
}

stark::core::NewtonsMethod::EvalResult stark::core::NewtonsMethod::_evaluate_E_grad_hess()
{
	this->_run_before_evaluation();
	this->logger->start_timing("evaluate_E_grad_hess");

	const int ndofs = this->global_energy->get_total_n_dofs();
	const int n_threads = this->settings->execution.n_threads;

	double E = 0.0;
	Eigen::VectorXd grad(ndofs);
	grad.setZero();
	auto element_hessians = this->compiled_global->evaluate_P__dP_du__local_d2P_du2(E, grad);

	if (this->settings->newton.project_to_PD) {
		element_hessians->project_to_PD_inplace__all(1e-12, /*use_mirroring=*/false);
	}
	auto hess = element_hessians->assemble_global(n_threads, ndofs);

	this->logger->stop_timing_add("evaluate_E_grad_hess");
	this->_run_after_evaluation();

	return { E, std::move(grad), hess };
}

stark::core::NewtonsMethod::EvalResult stark::core::NewtonsMethod::_evaluate_E_grad()
{
	this->_run_before_evaluation();
	this->logger->start_timing("evaluate_E_grad");

	const int ndofs = this->global_energy->get_total_n_dofs();

	double E = 0.0;
	Eigen::VectorXd grad(ndofs);
	grad.setZero();
	this->compiled_global->evaluate_P__dP_du(E, grad);

	this->logger->stop_timing_add("evaluate_E_grad");
	this->_run_after_evaluation();

	return { E, std::move(grad), nullptr };
}

stark::core::NewtonsMethod::EvalResult stark::core::NewtonsMethod::_evaluate_E()
{
	this->_run_before_evaluation();
	this->logger->start_timing("evaluate_E");

	double E = 0.0;
	this->compiled_global->evaluate_P(E);

	this->logger->stop_timing_add("evaluate_E");
	this->_run_after_evaluation();

	return { E, {}, nullptr };
}

Eigen::VectorXd stark::core::NewtonsMethod::_compute_residual(const Eigen::VectorXd& grad, double dt)
{
	if (this->settings->newton.residual.type == ResidualType::Force) {
		return grad.cwiseAbs() / dt;
	}
	else if (this->settings->newton.residual.type == ResidualType::Acceleration) {
		const std::vector<int32_t>& dofs_offsets = this->global_energy->get_dofs_offsets();
		this->residual = grad.cwiseAbs() / dt;
		int dof_count = 0;

		auto& inv_mass_application = this->callbacks->get_inv_mass_application();
		for (auto& pair : inv_mass_application) {
			const int dof = pair.first;
			auto& f = pair.second;

			const int begin = dofs_offsets[dof];
			const int end = dofs_offsets[dof + 1];
			dof_count += end - begin;
			f(this->residual.data() + begin, this->residual.data() + end);
		}

		if (dof_count != (int)this->residual.size()) {
			std::cout << "Stark error: NewtonsMethod::_compute_residual() found that not all dofs were used for ResidualType::Acceleration." << std::endl;
			exit(-1);
		}

		return this->residual;
	}
	else {
		std::cout << "Stark error: NewtonsMethod::_compute_residual() found an unknown residual type." << std::endl;
		exit(-1);
	}
}

double stark::core::NewtonsMethod::_compute_acceleration_correction(double du, double dt)
{
	return du / dt;
}

double stark::core::NewtonsMethod::_forcing_sequence(const Eigen::VectorXd& rhs)
{
	const double grad_norm = rhs.norm();
	return std::min(0.1, grad_norm * std::min(0.5, std::sqrt(grad_norm)));
}

typedef Eigen::SparseMatrix<double, EIGEN_DEFAULT_MATRIX_STORAGE_ORDER_OPTION, long long> longSparseMatrix;

bool stark::core::NewtonsMethod::_solve_linear_system(Eigen::VectorXd& du, const EvalResult& result, double dt)
{
	Eigen::VectorXd rhs = -1.0 * result.grad;

	switch (this->settings->newton.linear_system_solver)
	{
	case LinearSystemSolver::CG:
		{
			this->logger->start_timing("CG");
			const int n_threads = this->settings->execution.n_threads;
			const double cg_tol = this->_forcing_sequence(rhs);
			const int max_iterations = std::max(1000, (int)(this->settings->newton.cg_max_iterations_multiplier * rhs.size()));
			auto& lhs = *result.hess;
			this->du.resize(rhs.size());
			this->du.setZero();

			lhs.set_preconditioner(bsm::Preconditioner::BlockDiagonal);
			lhs.prepare_preconditioning(n_threads);
			static bsm::PCGContext pcg_ctx;
			bsm::PCGInfo info = bsm::solve_pcg(lhs, this->du.data(), rhs.data(), (int)rhs.size(), cg_tol, cg_tol, max_iterations, n_threads, false, pcg_ctx);
			this->cg_iterations_in_step += info.n_iterations;
			this->logger->stop_timing_add("CG");

			if (!info.converged) {
				this->console->add_error_msg(fmt::format("CG didn't converge to tolerance {:.2e} in max iterations {:d}.", cg_tol, max_iterations));
				return false;
			}
			return true;
		}
	case LinearSystemSolver::DirectLU:
		{
			this->logger->start_timing("directLU");
			std::vector<Eigen::Triplet<double>> triplets;
			result.hess->to_triplets(triplets);

			Eigen::SparseMatrix<double> s;
			s.resize(rhs.size(), rhs.size());
			s.setFromTriplets(triplets.begin(), triplets.end());
			s.makeCompressed();

			Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> lu;
			lu.analyzePattern(s);
			lu.factorize(s);
			du = lu.solve(rhs);
			this->logger->stop_timing_add("directLU");

			if (lu.info() != Eigen::ComputationInfo::Success) {
				this->console->add_error_msg("DirectLU couldn't find a solution.");
				return false;
			}
			return true;
		}
	case LinearSystemSolver::MKL_LU:
		{
#ifdef STARK_HAS_MKL
			this->logger->start_timing("MKL_LU");
			std::vector<Eigen::Triplet<double>> triplets;
			result.hess->to_triplets(triplets);

			longSparseMatrix s;
			s.resize(rhs.size(), rhs.size());
			s.setFromTriplets(triplets.begin(), triplets.end());
			s.makeCompressed();

			Eigen::PardisoLU<longSparseMatrix> lu(s);
			lu.analyzePattern(s);
			lu.factorize(s);
			du = lu.solve(rhs);
			this->logger->stop_timing_add("MKL_LU");

			if (lu.info() != Eigen::ComputationInfo::Success) {
				this->console->add_error_msg("MKL_LU couldn't find a solution.");
				return false;
			}
			return true;
#else
			this->console->add_error_msg("MKL_LU solver requested but MKL was not found at build time.");
			return false;
#endif
		}
	case LinearSystemSolver::MKL_LDLT:
		{
#ifdef STARK_HAS_MKL
			this->logger->start_timing("MKL_LDLT");
			std::vector<Eigen::Triplet<double>> triplets;
			result.hess->to_triplets(triplets);

			longSparseMatrix s;
			s.resize(rhs.size(), rhs.size());
			s.setFromTriplets(triplets.begin(), triplets.end());
			s.makeCompressed();

			Eigen::PardisoLDLT<longSparseMatrix> lu(s);
			lu.analyzePattern(s);
			lu.factorize(s);
			du = lu.solve(rhs);
			this->logger->stop_timing_add("MKL_LDLT");

			if (lu.info() != Eigen::ComputationInfo::Success) {
				this->console->add_error_msg("MKL_LDLT couldn't find a solution.");
				return false;
			}
			return true;
#else
			this->console->add_error_msg("MKL_LDLT solver requested but MKL was not found at build time.");
			return false;
#endif
		}
	default:
		std::cout << "Stark error: NewtonsMethod::_solve_linear_system() found an unknown linear system solver." << std::endl;
		exit(-1);
	}
}

double stark::core::NewtonsMethod::_inplace_max_step_in_search_direction(const Eigen::VectorXd& du)
{
	double step = this->callbacks->run_max_allowed_step();

	this->global_energy->get_dofs(this->u0.data());
	while (true) {
		if (step < 0.01) {
			return 0.0;
		}

		this->u1 = this->u0 + step * this->du;
		this->global_energy->set_dofs(this->u1.data());

		this->logger->start_timing("is_intermidiate_state_valid");
		const bool is_valid_state = this->callbacks->run_is_intermidiate_state_valid();
		this->logger->stop_timing_add("is_intermidiate_state_valid");

		if (is_valid_state) {
			break;
		}
		step *= 0.5;
	}
	this->console->print(fmt::format("max step = {:.2e} | ", step), ConsoleVerbosity::NewtonIterations);
	return step;
}

double stark::core::NewtonsMethod::_inplace_backtracking_line_search(const Eigen::VectorXd& du, double E0, double E, double step_valid_configuration, double du_dot_grad)
{
	double step = step_valid_configuration;
	const double suitable_backtracking_energy = E0 + 1e-4 * du_dot_grad;
	this->logger->start_timing("line_search");
	int line_search_it = 1;
	while (E > suitable_backtracking_energy) {
		this->step_line_search_count++;

		this->console->print(fmt::format("\n\t\t\t {:d}. step = {:.2e} | E/E_bt = {:.2e}", line_search_it, step, E / suitable_backtracking_energy), ConsoleVerbosity::NewtonIterations);

		if (line_search_it == this->settings->newton.max_line_search_iterations) {
			step = 0.0;
			break;
		}

		step *= this->settings->newton.line_search_multiplier;
		this->u1 = this->u0 + step * this->du;
		this->global_energy->set_dofs(this->u1.data());

		E = this->_evaluate_E().E;

		line_search_it++;
	}
	this->logger->stop_timing_add("line_search");

	if (step == 0 && this->settings->debug.line_search_output) {
		const std::string label = std::to_string(this->debug_output_counter);
		this->line_search_debug_logger.append_to_series(label, fmt::format("{:.6e}", E0));
		this->line_search_debug_logger.append_to_series(label, fmt::format("{:.6e}", suitable_backtracking_energy / E0));
		this->line_search_debug_logger.append_to_series(label, fmt::format("{:.6e}", this->du.norm()));
		for (double fstep = -1.0; fstep < 2.0; fstep += 0.01) {
			if (this->debug_output_counter == 0) {
				this->line_search_debug_logger.append_to_series("normalized_step_length", fmt::format("{:.6e}", fstep * step_valid_configuration));
			}
			this->u1 = this->u0 + fstep * step_valid_configuration * this->du;
			this->global_energy->set_dofs(this->u1.data());
			this->callbacks->run_before_energy_evaluation();
			double E_sample = 0.0;
			this->compiled_global->evaluate_P(E_sample);
			this->callbacks->run_after_energy_evaluation();
			this->line_search_debug_logger.append_to_series(label, fmt::format("{:.6e}", E_sample / E0));
		}
		this->line_search_debug_logger.save_to_disk();
		this->debug_output_counter++;
	}

	return step;
}
