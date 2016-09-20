//#pragma once
//#include "DeBoorKnotsGenerator.h"
//#include "Tridiagonal.h"
//
//namespace splineknots
//{
//	class EnhancedReducedDeboorKnotsGenerator final
//	{
//		InterpolativeMathFunction function_;
//		Tridiagonals tridagonals_;
//		bool is_parallel_;
//		DeBoorKnotsGenerator deboor_;
//	public:
//		EnhancedReducedDeboorKnotsGenerator(MathFunction function,
//			bool buffered = true);
//		EnhancedReducedDeboorKnotsGenerator(InterpolativeMathFunction function,
//			bool buffered = true);
//		KnotMatrix GenerateKnots(const SurfaceDimension& udimension,
//			const SurfaceDimension& vdimension,
//			double* calculation_time = nullptr);
//
//		void InParallel(bool in_parallel);
//
//		Tridiagonals& Tridagonals();
//
//		Tridiagonal& Tridiagonal();
//
//		struct PrecalculatedEnhancedReduced final
//		{
//			DeBoorKnotsGenerator::Precalculated deboor_precalculated_;
//			double twelve_div_h,
//				three_h_div_4;
//
//			PrecalculatedEnhancedReduced(const double h);
//		};
//
//
//		void InitializeKnots(const SurfaceDimension& udimension,
//			const SurfaceDimension& vdimension, KnotMatrix& values);
//
//		void Precalculate(const SurfaceDimension& udimension,
//			const SurfaceDimension& vdimension);
//
//		void InitializeBuffers(size_t u_count, size_t v_count);
//
//		template <typename RightSideSelector>
//		void RightSide(const RightSideSelector& right_side_variables,
//			const PrecalculatedEnhancedReduced& precalculated,
//			const double dfirst, const double dlast,
//			const int unknowns_count, KnotVector& rightside)
//		{
//			auto even = unknowns_count % 2 == 0;
//			auto tau = even ? 0 : 2;
//			auto eta = even ? -4 : 1;
//			auto upsilon = even ? unknowns_count : unknowns_count - 1;
//			auto equations_count = even ? unknowns_count / 2 - 1
//				: unknowns_count / 2;
//			auto three_div_h = precalculated.deboor_precalculated_.three_div_h;
//			auto twelve_div_h = three_div_h * 4;
//
//			rightside[0] = three_div_h * (right_side_variables(4)
//				- right_side_variables(0))
//				- twelve_div_h * (right_side_variables(3)
//					- right_side_variables(1)) - dfirst;
//
//			rightside[equations_count - 1] = three_div_h *
//				(right_side_variables(upsilon + tau) -
//					right_side_variables(upsilon - 2))
//				-
//				twelve_div_h *
//				(right_side_variables(upsilon + 1)
//					- right_side_variables(upsilon - 1)) - eta * dlast;
//
//			for (auto k = 2; k < equations_count; k++)
//			{
//				auto k2 = k * 2;
//				rightside[k - 1] =
//					three_div_h * (right_side_variables(2 * (k + 1))
//						- right_side_variables(2 * (k - 1))
//						- twelve_div_h * (right_side_variables(k2 + 1)
//							- right_side_variables(k2 - 1)));
//			}
//
//			// I do not know (yet) why but these must be half of values
//			// designed by L. Mino.
//			// Therefore this loop shouldn't be here
//			for (int i = 1; i < equations_count - 1; i++)
//			{
//				rightside[i] *= 0.5;
//			}
//		}
//
//		void RightSideCross(const KnotMatrix& knots, const int i,
//			const double dfirst, const double dlast,
//			const int unknowns_count, KnotVector& rightside);
//
//		void FillXDerivations(KnotMatrix& values);
//
//		void FillXYDerivations(KnotMatrix& values);
//
//		void FillYDerivations(KnotMatrix& values);
//
//		void FillYXDerivations(KnotMatrix& values);
//
//		void FillXDerivations(const int column_index, KnotMatrix& values);
//
//		void FillXYDerivations(const int column_index, KnotMatrix& values);
//
//		void FillYDerivations(const int row_index, KnotMatrix& values);
//
//		void FillYXDerivations(const int row_index, KnotMatrix& values);
//
//		template <typename RightSideSelector, typename UnknownsSetter>
//		void SolveTridiagonal(const RightSideSelector& selector,
//			const PrecalculatedEnhancedReduced& precalculated, const double dfirst,
//			const double dlast, const int unknowns_count,
//			UnknownsSetter& unknowns_setter)
//		{
//			auto& tridiagonal = Tridiagonal();
//			auto& rightside = Tridiagonal().RightSideBuffer();
//			RightSide(selector, precalculated, dfirst, dlast, unknowns_count, rightside);
//			auto& result = tridiagonal.Solve(unknowns_count);
//			for (size_t k = 0; k < unknowns_count / 2 - 1; k++)
//			{
//				unknowns_setter(2 * (k + 1), result[k]);
//			}
//		}
//
//	private:
//		PrecalculatedEnhancedReduced precalculated_hx_;
//		PrecalculatedEnhancedReduced precalculated_hy_;
//
//	};
//}
