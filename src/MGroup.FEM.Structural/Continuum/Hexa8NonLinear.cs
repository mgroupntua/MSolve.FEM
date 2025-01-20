namespace MGroup.FEM.Structural.Continuum
{
	using System;
	using System.Collections.Generic;
	using System.Linq;

	using MGroup.Constitutive.Structural;
	using MGroup.Constitutive.Structural.Continuum;
	using MGroup.LinearAlgebra.Implementations;
	using MGroup.LinearAlgebra.Matrices;
	using MGroup.MSolve.Constitutive;
	using MGroup.MSolve.DataStructures;
	using MGroup.MSolve.Discretization;
	using MGroup.MSolve.Discretization.BoundaryConditions;
	using MGroup.MSolve.Discretization.Dofs;
	using MGroup.MSolve.Discretization.Embedding;
	using MGroup.MSolve.Discretization.Entities;
	using MGroup.MSolve.Numerics.Integration.Quadratures;
	using MGroup.MSolve.Numerics.Interpolation;
	using MGroup.MSolve.Numerics.Interpolation.Jacobians;

	/// <summary>
	/// Continuum finite Element for 3d problems with material and geometric nonlinearities
	/// Authors: Gerasimos Sotiropoulos
	/// </summary>
	public class Hexa8NonLinear : IStructuralElementType, IEmbeddedHostElement
	{
		protected readonly static IDofType[] nodalDOFTypes = new IDofType[] { StructuralDof.TranslationX, StructuralDof.TranslationY, StructuralDof.TranslationZ };
		protected readonly static IDofType[][] dofTypes = new IDofType[][] { nodalDOFTypes, nodalDOFTypes, nodalDOFTypes,
			nodalDOFTypes, nodalDOFTypes, nodalDOFTypes, nodalDOFTypes, nodalDOFTypes };
		protected readonly IContinuumMaterial3D[] materialsAtGaussPoints;
		protected IElementDofEnumerator dofEnumerator = new GenericDofEnumerator();

		private readonly int nGaussPoints;
		private bool isInitialized = false;

		private double[][] initialCoordinates; //not defined by user. 8 arrays of 3 elements
		private double[][] totalDisplacements;
		private double[] integrationCoeffs;

		private double[][] strainsVec;
		private double[][] strainsVec_last_converged;
		private double[][] lastStresses;

		protected Hexa8NonLinear()
		{
		}

		public Hexa8NonLinear(IReadOnlyList<INode> nodes, IContinuumMaterial3D material, IQuadrature3D quadratureForStiffness)
		{
			this.nGaussPoints = quadratureForStiffness.IntegrationPoints.Count;
			this.QuadratureForStiffness = quadratureForStiffness;
			this.Interpolation = InterpolationHexa8.UniqueInstance;
			Nodes = nodes;
			lastStresses = new double[nGaussPoints][];
			materialsAtGaussPoints = new IContinuumMaterial3D[nGaussPoints];
			for (int i = 0; i < nGaussPoints; i++)
			{
				materialsAtGaussPoints[i] = (IContinuumMaterial3D)material.Clone();
				lastStresses[i] = new double[6];
			}
		}

		public InterpolationHexa8 Interpolation { get; }
		public IQuadrature3D QuadratureForStiffness { get; }

		public int ID { get; set; }
		public int SubdomainID { get; set; }
		public IReadOnlyList<INode> Nodes { get; }
		public CellType CellType { get; } = CellType.Hexa8;

		public ElementDimensions ElementDimensions => ElementDimensions.ThreeD;

		public IElementDofEnumerator DofEnumerator
		{
			get { return dofEnumerator; }
			set { dofEnumerator = value; }
		}

		//public bool ConstitutiveLawModified
		//{
		//	get
		//	{
		//		foreach (IContinuumMaterial3D material in materialsAtGaussPoints)
		//			if (material.IsCurrentStateDifferent()) return true;
		//		return false;
		//	}
		//}

		private Matrix[] Getbl13DeformationMatrices(IReadOnlyList<Matrix> shapeFunctionNaturalDerivatives)
		{
			Matrix[] bl13Matrices;
			bl13Matrices = new Matrix[nGaussPoints];
			for (int npoint = 0; npoint < nGaussPoints; npoint++)
			{
				bl13Matrices[npoint] = Matrix.CreateZero(9, 24);
				for (int m = 0; m < 8; m++)
				{
					for (int n = 0; n < 3; n++)
					{
						bl13Matrices[npoint][n, 3 * m + 0] = shapeFunctionNaturalDerivatives[npoint][m,n];
						bl13Matrices[npoint][n + 3, 3 * m + 1] = shapeFunctionNaturalDerivatives[npoint][m,n];
						bl13Matrices[npoint][n + 6, 3 * m + 2] = shapeFunctionNaturalDerivatives[npoint][m,n];
					}
				}
			}
			return bl13Matrices;
		}

		private Matrix[] Getbl11aDeformationMatrices(Matrix[] jacobianInverse)
		{
			Matrix[] bl11aMatrices = new Matrix[nGaussPoints];
			for (int gpoint = 0; gpoint < nGaussPoints; gpoint++)
			{
				bl11aMatrices[gpoint] = Matrix.CreateZero(6, 9);
				for (int m = 0; m < 3; m++) // calculate first three data lines of the matrix
				{
					for (int n = 0; n < 3; n++)
					{
						bl11aMatrices[gpoint][m, 3 * m + n] = jacobianInverse[gpoint][m, n];
					}
				}
				for (int n = 0; n < 3; n++)
				{
					bl11aMatrices[gpoint][3, n] = jacobianInverse[gpoint][1, n]; // calculate 4th data line
					bl11aMatrices[gpoint][3, 3 + n] = jacobianInverse[gpoint][0, n];
					bl11aMatrices[gpoint][4, 3 + n] = jacobianInverse[gpoint][2, n]; // calculate 5th data line
					bl11aMatrices[gpoint][4, 6 + n] = jacobianInverse[gpoint][1, n];
					bl11aMatrices[gpoint][5, 0 + n] = jacobianInverse[gpoint][2, n]; // calculate 6th data line
					bl11aMatrices[gpoint][5, 6 + n] = jacobianInverse[gpoint][0, n];
				}
			}

			return bl11aMatrices;
		}

		private Matrix[] GetBL12DeformationMatrices(Matrix[] jacobianInverse)
		{
			Matrix[] bl12Marices = new Matrix[nGaussPoints];
			for (int gpoint = 0; gpoint < nGaussPoints; gpoint++)
			{
				bl12Marices[gpoint] = Matrix.CreateZero(9, 9);
				for (int m = 0; m < 3; m++) // calculate first three data lines of the matrix
				{
					for (int n = 0; n < 3; n++)
					{
						bl12Marices[gpoint][m, 3 * m + n] = jacobianInverse[gpoint][0, n];
					}
				}
				for (int m = 0; m < 3; m++) // calculate  data lines 4:6
				{
					for (int n = 0; n < 3; n++)
					{
						bl12Marices[gpoint][3 + m, 3 * m + n] = jacobianInverse[gpoint][1, n];
					}
				}
				for (int m = 0; m < 3; m++) // calculate  data lines 7:8
				{
					for (int n = 0; n < 3; n++)
					{
						bl12Marices[gpoint][6 + m, 3 * m + n] = jacobianInverse[gpoint][2, n];
					}
				}

			}

			return bl12Marices;
		}

		private Matrix[] Getbl01MDeformationMatrices(Matrix[] jacobianInverse)
		{
			Matrix[] bl01Matrices = new Matrix[nGaussPoints];
			for (int gpoint = 0; gpoint < nGaussPoints; gpoint++)
			{
				bl01Matrices[gpoint] = Matrix.CreateZero(6, 9);
				for (int m = 0; m < 3; m++) // calculate first three data lines of the matrix
				{
					for (int n = 0; n < 3; n++)
					{
						bl01Matrices[gpoint][m, 3 * m + n] = jacobianInverse[gpoint][m, n];
					}
				}
				for (int n = 0; n < 3; n++)
				{
					bl01Matrices[gpoint][3, n] = jacobianInverse[gpoint][1, n]; // calculate 4th data line
					bl01Matrices[gpoint][3, 3 + n] = jacobianInverse[gpoint][0, n];
					bl01Matrices[gpoint][4, 3 + n] = jacobianInverse[gpoint][2, n]; // calculate 5th data line
					bl01Matrices[gpoint][4, 6 + n] = jacobianInverse[gpoint][1, n];
					bl01Matrices[gpoint][5, 0 + n] = jacobianInverse[gpoint][2, n]; // calculate 6th data line
					bl01Matrices[gpoint][5, 6 + n] = jacobianInverse[gpoint][0, n];
				}
			}
			return bl01Matrices;
		}

		private Matrix[] GetAuxilliaryDeformationbnl1Matrices(Matrix[] jacobianInverse)
		{
			Matrix[] bnl1Matrices = new Matrix[nGaussPoints];
			for (int gpoint = 0; gpoint < nGaussPoints; gpoint++)
			{
				bnl1Matrices[gpoint] = Matrix.CreateZero(9, 9);
				for (int m = 0; m < 3; m++)
				{
					for (int n = 0; n < 3; n++)
					{
						for (int p = 0; p < 3; p++)
						{
							bnl1Matrices[gpoint][3 * m + n, 3 * m + p] = jacobianInverse[gpoint][n, p];
						}
					}
				}
			}
			return bnl1Matrices;
		}

		private void CalculateInitialConfigurationData()
		{
			IReadOnlyList<Matrix> shapeFunctionNaturalDerivatives;
			shapeFunctionNaturalDerivatives = Interpolation.EvaluateNaturalGradientsAtGaussPoints(QuadratureForStiffness);

			Matrix[] bl13Matrices;
			bl13Matrices = Getbl13DeformationMatrices(shapeFunctionNaturalDerivatives);

			Matrix[] bnl1Matrices;

			initialCoordinates = new double[8][];
			totalDisplacements = new double[8][];

			(Matrix[] jacobianInverse, double[] jacobianDeterminants) = JacobianHexa8Reverse.GetJ_0invHexaAndjacobianDeterminants(
				shapeFunctionNaturalDerivatives, Nodes, nGaussPoints);

			integrationCoeffs = new double[nGaussPoints];

			bnl1Matrices = GetAuxilliaryDeformationbnl1Matrices(jacobianInverse);

			for (int j = 0; j < 8; j++)
			{
				initialCoordinates[j] = new double[] { Nodes[j].X, Nodes[j].Y, Nodes[j].Z, };
			}

			for (int gpoint = 0; gpoint < nGaussPoints; gpoint++)
			{
				integrationCoeffs[gpoint] = jacobianDeterminants[gpoint] * QuadratureForStiffness.IntegrationPoints[gpoint].Weight;

			}

			totalDisplacements = new double[8][];
			strainsVec = new double[nGaussPoints][];
			strainsVec_last_converged = new double[nGaussPoints][];
			for (int gpoint = 0; gpoint < nGaussPoints; gpoint++)
			{
				strainsVec[gpoint] = new double[6];
				strainsVec_last_converged[gpoint] = new double[6];
			}
			for (int k = 0; k < 8; k++)
			{
				totalDisplacements[k] = new double[3];
			}
			isInitialized = true;

		}

		private void UpdateCoordinateData(double[] localdisplacements, out double[][] deformedCoordinates)
		{
			deformedCoordinates = new double[8][];
			for (int j = 0; j < 8; j++)
			{
				deformedCoordinates[j] = new double[3];
				for (int k = 0; k < 3; k++)
				{
					totalDisplacements[j][k] = localdisplacements[3 * j + k];
					deformedCoordinates[j][k] = initialCoordinates[j][k] + totalDisplacements[j][k];
				}
			}
		}

		private void CalculateStrains(double[] localdisplacements, double[][] deformedCoordinates)
		{
			IReadOnlyList<Matrix> shapeFunctionNaturalDerivatives;
			shapeFunctionNaturalDerivatives = Interpolation.EvaluateNaturalGradientsAtGaussPoints(QuadratureForStiffness);
			(Matrix[] jacobianInverse, double[] jacobianDeterminants) = JacobianHexa8Reverse.GetJ_0invHexaAndjacobianDeterminants(
				shapeFunctionNaturalDerivatives, Nodes, nGaussPoints);
			//TODO: possibility of caching shapeFunctionNaturalDerivatives or J_0inv

			Matrix[] deformationGradientsTransposed = new Matrix[nGaussPoints];
			Matrix[] GL = new Matrix[nGaussPoints];
			for (int npoint = 0; npoint < nGaussPoints; npoint++)
			{
				deformationGradientsTransposed[npoint] = Matrix.CreateZero(3, 3);
				GL[npoint] = Matrix.CreateZero(3, 3);
			}

			Matrix[] jacobiansDeformedMatrices = JacobianHexa8Reverse.Get_jacobiansDeformedMatrices(nGaussPoints, deformedCoordinates, shapeFunctionNaturalDerivatives);

			//
			for (int npoint = 0; npoint < nGaussPoints; npoint++)
			{
				//
				deformationGradientsTransposed[npoint] = jacobianInverse[npoint] * jacobiansDeformedMatrices[npoint];

				//
				GL[npoint] = deformationGradientsTransposed[npoint] * deformationGradientsTransposed[npoint].Transpose();
				for (int m = 0; m < 3; m++)
				{
					GL[npoint][m, m] += -1;
				}
				GL[npoint].ScaleIntoThis(0.5);

				//
				for (int m = 0; m < 3; m++)
				{
					strainsVec[npoint][m] = GL[npoint][m, m];
				}
				strainsVec[npoint][3] = 2 * GL[npoint][0, 1];
				strainsVec[npoint][4] = 2 * GL[npoint][1, 2];
				strainsVec[npoint][5] = 2 * GL[npoint][2, 0];
			}

		}

		private double[] UpdateResponseIntegral()
		{
			//TODO: the gauss point loop should be the outer one

			// Matrices that are not currently cached are calculated here.
			Matrix ll2 = Matrix.CreateZero(8, 3);
			for (int m = 0; m < 8; m++)
			{
				for (int n = 0; n < 3; n++)
				{
					ll2[m, n] = totalDisplacements[m][n];
				}
			}
			IReadOnlyList<Matrix> shapeFunctionNaturalDerivatives;
			shapeFunctionNaturalDerivatives = Interpolation.EvaluateNaturalGradientsAtGaussPoints(QuadratureForStiffness);
			(Matrix[] jacobianInverse, double[] jacobianDeterminants) = JacobianHexa8Reverse.GetJ_0invHexaAndjacobianDeterminants(
				shapeFunctionNaturalDerivatives, Nodes, nGaussPoints);
			Matrix[] bl13Matrices;
			bl13Matrices = Getbl13DeformationMatrices(shapeFunctionNaturalDerivatives);
			Matrix[] bl11aMatrices; // dimension number of gpoints
			Matrix[] bl12Marices;
			Matrix[] bl01Matrices;
			bl11aMatrices = Getbl11aDeformationMatrices(jacobianInverse);
			bl12Marices = GetBL12DeformationMatrices(jacobianInverse);
			bl01Matrices = Getbl01MDeformationMatrices(jacobianInverse);

			//INITIALIZATION of MAtrixes that are currently not cached
			double[][] integrCoeffsTimesStresses = new double[nGaussPoints][];
			Matrix[] BL = new Matrix[nGaussPoints];
			for (int gpoint = 0; gpoint < nGaussPoints; gpoint++)
			{
				integrCoeffsTimesStresses[gpoint] = new double[6];
				BL[gpoint] = Matrix.CreateZero(6, 24);
			}

			double[][] forces = new double[nGaussPoints + 1][];
			for (int npoint = 0; npoint < nGaussPoints + 1; npoint++)
			{
				forces[npoint] = new double[24];
			}

			Matrix[] BL11 = new Matrix[nGaussPoints];
			Matrix[] bL1112Plus01Matrices = new Matrix[nGaussPoints];
			for (int npoint = 0; npoint < nGaussPoints; npoint++)
			{
				BL11[npoint] = Matrix.CreateZero(6, 9);
				bL1112Plus01Matrices[npoint] = Matrix.CreateZero(6, 9);
			}

			for (int npoint = 0; npoint < nGaussPoints; npoint++)
			{

				integrCoeffsTimesStresses[npoint] = lastStresses[npoint].Scale(integrationCoeffs[npoint]);

				//
				Matrix lcyrcumflex = Matrix.CreateZero(3, 3);
				lcyrcumflex = shapeFunctionNaturalDerivatives[npoint].Transpose() * ll2;

				for (int m = 0; m < 6; m++)
				{
					for (int n = 0; n < 3; n++)
					{
						for (int p = 0; p < 3; p++)
						{
							BL11[npoint][m, n] += bl11aMatrices[npoint][m, p] * lcyrcumflex[p, n];
							BL11[npoint][m, 3 + n] += bl11aMatrices[npoint][m, 3 + p] * lcyrcumflex[p, n];
							BL11[npoint][m, 6 + n] += bl11aMatrices[npoint][m, 6 + p] * lcyrcumflex[p, n];
						}
					}
				}

				//
				bL1112Plus01Matrices[npoint] = BL11[npoint] * bl12Marices[npoint];
				bL1112Plus01Matrices[npoint].AddIntoThis(bl01Matrices[npoint]);

				// 
				BL[npoint] = bL1112Plus01Matrices[npoint] * bl13Matrices[npoint];

				//              
				forces[npoint] = BL[npoint].Multiply(integrCoeffsTimesStresses[npoint], true);
			}

			for (int npoint = 0; npoint < nGaussPoints; npoint++)
			{
				forces[nGaussPoints].AddIntoThis(forces[npoint]);
			}

			return forces[nGaussPoints];
		}

		private Matrix UpdateKmatrices()
		{
			Matrix elementStiffnessMatrix = Matrix.CreateZero(24, 24);

			// initialization of matrices that are not cached currently
			double[][] integrCoeffsTimesStresses = new double[nGaussPoints][];
			Matrix[] BL = new Matrix[nGaussPoints];
			for (int gpoint = 0; gpoint < nGaussPoints; gpoint++)
			{
				integrCoeffsTimesStresses[gpoint] = new double[6];
				BL[gpoint] = Matrix.CreateZero(6, 24);

			}

			Matrix ll2 = Matrix.CreateZero(8, 3);
			for (int m = 0; m < 8; m++)
			{
				for (int n = 0; n < 3; n++)
				{
					ll2[m, n] = totalDisplacements[m][n];
				}
			}

			IReadOnlyList<Matrix> shapeFunctionNaturalDerivatives;
			shapeFunctionNaturalDerivatives = Interpolation.EvaluateNaturalGradientsAtGaussPoints(QuadratureForStiffness);
			(Matrix[] jacobianInverse, double[] jacobianDeterminants) = JacobianHexa8Reverse.GetJ_0invHexaAndjacobianDeterminants(
				shapeFunctionNaturalDerivatives, Nodes, nGaussPoints);
			Matrix[] bl13Matrices;
			bl13Matrices = Getbl13DeformationMatrices(shapeFunctionNaturalDerivatives);
			Matrix[] bl11aMatrices; // dimension: gpoints
			Matrix[] bl12Marices;
			Matrix[] bl01Matrices;
			bl11aMatrices = Getbl11aDeformationMatrices(jacobianInverse);
			bl12Marices = GetBL12DeformationMatrices(jacobianInverse);
			bl01Matrices = Getbl01MDeformationMatrices(jacobianInverse);

			Matrix[] BL11 = new Matrix[nGaussPoints];
			Matrix[] bL1112Plus01Matrices = new Matrix[nGaussPoints];
			for (int npoint = 0; npoint < nGaussPoints; npoint++)
			{
				BL11[npoint] = Matrix.CreateZero(6, 9);
				bL1112Plus01Matrices[npoint] = Matrix.CreateZero(6, 9); //TODO this may be unnescessary
			}

			var s = MatrixSymmetry.Symmetric;
			for (int npoint = 0; npoint < nGaussPoints; npoint++)
			{

				// 
				integrCoeffsTimesStresses[npoint] = lastStresses[npoint].Scale(integrationCoeffs[npoint]);

				//
				Matrix lcyrcumflex = Matrix.CreateZero(3, 3);
				lcyrcumflex = shapeFunctionNaturalDerivatives[npoint].Transpose() * ll2;

				for (int m = 0; m < 6; m++)
				{
					for (int n = 0; n < 3; n++)
					{
						for (int p = 0; p < 3; p++)
						{
							BL11[npoint][m, n] += bl11aMatrices[npoint][m, p] * lcyrcumflex[p, n];
							BL11[npoint][m, 3 + n] += bl11aMatrices[npoint][m, 3 + p] * lcyrcumflex[p, n];
							BL11[npoint][m, 6 + n] += bl11aMatrices[npoint][m, 6 + p] * lcyrcumflex[p, n];
						}
					}
				}

				// 
				bL1112Plus01Matrices[npoint] = BL11[npoint] * bl12Marices[npoint];
				bL1112Plus01Matrices[npoint].AddIntoThis(bl01Matrices[npoint]);

				//
				BL[npoint] = bL1112Plus01Matrices[npoint] * bl13Matrices[npoint];

			}
			// TODO: BL and above calculations can cached from calculate forces method

			Matrix[] bnl1Matrices;
			Matrix[] bnlMatrices;
			bnl1Matrices = GetAuxilliaryDeformationbnl1Matrices(jacobianInverse);
			bnlMatrices = new Matrix[nGaussPoints];
			for (int gpoint = 0; gpoint < nGaussPoints; gpoint++)
			{
				bnlMatrices[gpoint] = Matrix.CreateZero(9, 24); //todo this may be unnescessary
				bnlMatrices[gpoint] = bnl1Matrices[gpoint] * bl13Matrices[gpoint];
			}

			Matrix[] integrCoeff_Spk = new Matrix[nGaussPoints];
			for (int npoint = 0; npoint < nGaussPoints; npoint++)
			{
				integrCoeff_Spk[npoint] = Matrix.CreateZero(3, 3);
			}

			Matrix[] kl_ = new Matrix[nGaussPoints + 1];
			Matrix[] knl_ = new Matrix[nGaussPoints + 1];
			for (int npoint = 0; npoint < nGaussPoints + 1; npoint++)
			{
				kl_[npoint] = Matrix.CreateZero(24, 24);
				knl_[npoint] = Matrix.CreateZero(24, 24);
			}

			for (int npoint = 0; npoint < nGaussPoints; npoint++)
			{
				Matrix integrCoeff_SPK_epi_bnlMatrices = Matrix.CreateZero(9, 24); //TODO
				Matrix integrCoeff_cons_disp = Matrix.CreateZero(6, 6); //TODO
				Matrix integrCoeff_cons_disp_epi_BL = Matrix.CreateZero(6, 24);//TODO

				//
				integrCoeff_Spk[npoint][0, 0] = integrCoeffsTimesStresses[npoint][0];
				integrCoeff_Spk[npoint][0, 1] = integrCoeffsTimesStresses[npoint][3];
				integrCoeff_Spk[npoint][0, 2] = integrCoeffsTimesStresses[npoint][5];
				integrCoeff_Spk[npoint][1, 0] = integrCoeffsTimesStresses[npoint][3];
				integrCoeff_Spk[npoint][1, 1] = integrCoeffsTimesStresses[npoint][1];
				integrCoeff_Spk[npoint][1, 2] = integrCoeffsTimesStresses[npoint][4];
				integrCoeff_Spk[npoint][2, 0] = integrCoeffsTimesStresses[npoint][5];
				integrCoeff_Spk[npoint][2, 1] = integrCoeffsTimesStresses[npoint][4];
				integrCoeff_Spk[npoint][2, 2] = integrCoeffsTimesStresses[npoint][2];

				//
				IMatrixView consDisp = materialsAtGaussPoints[npoint].ConstitutiveMatrix;
				s = s == MatrixSymmetry.Symmetric ? consDisp.MatrixSymmetry : s;

				for (int m = 0; m < 6; m++)
				{
					for (int n = 0; n < 6; n++)
					{
						integrCoeff_cons_disp[m, n] = integrationCoeffs[npoint] * consDisp[m, n];
					}
				}

				//
				integrCoeff_cons_disp_epi_BL = integrCoeff_cons_disp * BL[npoint];

				//
				kl_[npoint] = BL[npoint].Transpose() * integrCoeff_cons_disp_epi_BL;

				//
				for (int m = 0; m < 3; m++) // 3x24 dimensions
				{
					for (int n = 0; n < 24; n++)
					{
						for (int p = 0; p < 3; p++)
						{
							integrCoeff_SPK_epi_bnlMatrices[m, n] += integrCoeff_Spk[npoint][m, p] * bnlMatrices[npoint][p, n];
							integrCoeff_SPK_epi_bnlMatrices[3 + m, n] += integrCoeff_Spk[npoint][m, p] * bnlMatrices[npoint][3 + p, n];
							integrCoeff_SPK_epi_bnlMatrices[6 + m, n] += integrCoeff_Spk[npoint][m, p] * bnlMatrices[npoint][6 + p, n];
						}
					}
				}

				//
				knl_[npoint] = bnlMatrices[npoint].Transpose() * integrCoeff_SPK_epi_bnlMatrices;
			}

			// Add contributions of each gp on the total element stiffness matrix elementStiffnessMatrix            
			for (int npoint = 0; npoint < nGaussPoints; npoint++)
			{
				for (int m = 0; m < 24; m++)
				{
					for (int n = 0; n < 24; n++)
					{
						kl_[nGaussPoints][m, n] += kl_[npoint][m, n];
						knl_[nGaussPoints][m, n] += knl_[npoint][m, n];
					}
				}
			}
			for (int m = 0; m < 24; m++)
			{
				for (int n = 0; n < 24; n++)
				{
					elementStiffnessMatrix[m, n] = kl_[nGaussPoints][m, n] + knl_[nGaussPoints][m, n];
				}
			}

			elementStiffnessMatrix.MatrixSymmetry = s;
			return elementStiffnessMatrix;
		}

		public Tuple<double[], double[]> CalculateResponse(double[] localTotalDisplacements)
		{
			this.UpdateCoordinateData(localTotalDisplacements, out double[][] deformedCoordinates);
			CalculateStrains(localTotalDisplacements, deformedCoordinates);
			double[] strainsVec_strain_minus_last_converged_value = new double[6];
			for (int npoint = 0; npoint < materialsAtGaussPoints.Length; npoint++)
			{
				strainsVec_strain_minus_last_converged_value = new double[6]
				{
					strainsVec[npoint][0]- strainsVec_last_converged[npoint][0],
					strainsVec[npoint][1] - strainsVec_last_converged[npoint][1],
					strainsVec[npoint][2] - strainsVec_last_converged[npoint][2],
					strainsVec[npoint][3]- strainsVec_last_converged[npoint][3],
					strainsVec[npoint][4]- strainsVec_last_converged[npoint][4],
					strainsVec[npoint][5]- strainsVec_last_converged[npoint][5]
				};
				lastStresses[npoint] = materialsAtGaussPoints[npoint].UpdateConstitutiveMatrixAndEvaluateResponse(strainsVec_strain_minus_last_converged_value);
				//To update with total strain simplY = materialsAtGaussPoints[npoint].UpdateMaterial(strainsVec[npoint]);
			}
			return new Tuple<double[], double[]>(strainsVec_strain_minus_last_converged_value, lastStresses[materialsAtGaussPoints.Length - 1]);
			//TODO return data with total strains data would be:
			//return new Tuple<double[], double[]>(strainsVec[materialsAtGaussPoints.Length - 1], materialsAtGaussPoints[materialsAtGaussPoints.Length - 1].Stresses);
			//TODO: why return only the strain- stress of the gausspoint that is last on the array, Where is it needed?
		}

		public double[] CalculateResponseIntegral()
			=> this.UpdateResponseIntegral();

		public double[] CalculateResponseIntegralForLogging(double[] localDisplacements)
			=> CalculateResponseIntegral();

		public virtual IMatrix StiffnessMatrix()
		{
			if (!isInitialized)
			{
				this.CalculateInitialConfigurationData();
				var localTotalDisplacements = new double[24];
				this.UpdateCoordinateData(localTotalDisplacements, out double[][] deformedCoordinates);
				this.CalculateStrains(localTotalDisplacements, deformedCoordinates);
			}
			Matrix elementStiffness = this.UpdateKmatrices();
			//It doesn't implement Iembedded to return dof.Enumerator.GetTransformedMatrix
			return elementStiffness;
		}
		public IMatrix PhysicsMatrix()
		{
			return StiffnessMatrix();
		}
		//public void ResetConstitutiveLawModified()
		//{
		//	foreach (IContinuumMaterial3D material in materialsAtGaussPoints) material.ResetModified();
		//}

		public void ClearConstitutiveLawState()
		{
			//TODO: the next throws an exception. Investigate. Possible changes in Analyzers may be the cause.
			//foreach (IContinuumMaterial3D m in materialsAtGaussPoints) m.ClearState();
		}

		public void SaveConstitutiveLawState(IHaveState externalState)
		{
			for (int npoint = 0; npoint < materialsAtGaussPoints.Length; npoint++)
			{
				for (int i1 = 0; i1 < 6; i1++)
				{ strainsVec_last_converged[npoint][i1] = strainsVec[npoint][i1]; }
			}

			foreach (IContinuumMaterial3D m in materialsAtGaussPoints) m.CreateState();
			
			if (externalState != null && (externalState is IHaveStateWithValues))
			{
				var s = (IHaveStateWithValues)externalState;
				if (s.StateValues.ContainsKey(TransientLiterals.TIME))
				{
					var time = s.StateValues[TransientLiterals.TIME];
					foreach (var m in materialsAtGaussPoints.Where(x => x is ITransientConstitutiveLaw).Select(x => (ITransientConstitutiveLaw)x))
					{
						m.SetCurrentTime(time);
					}
				}

			}
		}

		//public void ClearMaterialStresses()
		//{
		//	foreach (IContinuumMaterial3D m in materialsAtGaussPoints) m.ClearStresses();
		//}

		public virtual IReadOnlyList<IReadOnlyList<IDofType>> GetElementDofTypes() => dofTypes;

		#region not implemented
		//public double[] CalculateAccelerationResponseIntegral(IElementType element, IList<MassAccelerationLoad> loads)
		//{
		//	throw new NotImplementedException();
		//}

		public virtual IMatrix MassMatrix()
		{
			throw new NotImplementedException();
		}

		public virtual IMatrix DampingMatrix()
		{
			throw new NotImplementedException();
		}
		#endregion


		#region IEmbeddedHostElement
		protected double[,] GetCoordinatesTranspose(IElementType element)
		{
			double[,] nodeCoordinatesXYZ = new double[3, dofTypes.Length];
			for (int i = 0; i < dofTypes.Length; i++)
			{
				nodeCoordinatesXYZ[0, i] = element.Nodes[i].X;
				nodeCoordinatesXYZ[1, i] = element.Nodes[i].Y;
				nodeCoordinatesXYZ[2, i] = element.Nodes[i].Z;
			}
			return nodeCoordinatesXYZ;
		}

		//TODO: This should be handled by InterpolationHexa8Reverse
		private double[] CalcH8Shape(double fXi, double fEta, double fZeta)
		{
			const double fSqC125 = 0.5;
			double auxilliaryfXiP = (1.0 + fXi) * fSqC125;
			double auxilliaryfEtaP = (1.0 + fEta) * fSqC125;
			double auxilliaryfZetaP = (1.0 + fZeta) * fSqC125;
			double auxilliaryfXiM = (1.0 - fXi) * fSqC125;
			double auxilliaryfEtaM = (1.0 - fEta) * fSqC125;
			double auxilliaryfZetaM = (1.0 - fZeta) * fSqC125;

			double[] auxH8ShapeFunctiondata = new double[8]; // Warning: shape function data not in hexa8fixed order.

			auxH8ShapeFunctiondata[6] = auxilliaryfXiM * auxilliaryfEtaM * auxilliaryfZetaM;
			auxH8ShapeFunctiondata[7] = auxilliaryfXiP * auxilliaryfEtaM * auxilliaryfZetaM;
			auxH8ShapeFunctiondata[4] = auxilliaryfXiP * auxilliaryfEtaP * auxilliaryfZetaM;
			auxH8ShapeFunctiondata[5] = auxilliaryfXiM * auxilliaryfEtaP * auxilliaryfZetaM;
			auxH8ShapeFunctiondata[2] = auxilliaryfXiM * auxilliaryfEtaM * auxilliaryfZetaP;
			auxH8ShapeFunctiondata[3] = auxilliaryfXiP * auxilliaryfEtaM * auxilliaryfZetaP;
			auxH8ShapeFunctiondata[0] = auxilliaryfXiP * auxilliaryfEtaP * auxilliaryfZetaP;
			auxH8ShapeFunctiondata[1] = auxilliaryfXiM * auxilliaryfEtaP * auxilliaryfZetaP;

			return auxH8ShapeFunctiondata;
		}

		//TODO: This should be handled by InterpolationHexa8Reverse
		private double[] CalcH8NablaShape(double fXi, double fEta, double fZeta)
		{
			const double fSq125 = 0.35355339059327376220042218105242;
			double auxilliaryfXiP = (1.0 + fXi) * fSq125;
			double auxilliaryfEtaP = (1.0 + fEta) * fSq125;
			double auxilliaryfZetaP = (1.0 + fZeta) * fSq125;
			double auxilliaryfXiM = (1.0 - fXi) * fSq125;
			double auxilliaryfEtaM = (1.0 - fEta) * fSq125;
			double auxilliaryfZetaM = (1.0 - fZeta) * fSq125;

			double[] auxilliaryfaDS = new double[24];
			auxilliaryfaDS[6] = -auxilliaryfEtaM * auxilliaryfZetaM;
			auxilliaryfaDS[4] = auxilliaryfEtaP * auxilliaryfZetaM;
			auxilliaryfaDS[2] = -auxilliaryfEtaM * auxilliaryfZetaP;
			auxilliaryfaDS[0] = auxilliaryfEtaP * auxilliaryfZetaP;
			auxilliaryfaDS[7] = -auxilliaryfaDS[6];
			auxilliaryfaDS[5] = -auxilliaryfaDS[4];
			auxilliaryfaDS[3] = -auxilliaryfaDS[2];
			auxilliaryfaDS[1] = -auxilliaryfaDS[0];


			auxilliaryfaDS[14] = -auxilliaryfXiM * auxilliaryfZetaM;
			auxilliaryfaDS[15] = -auxilliaryfXiP * auxilliaryfZetaM;
			auxilliaryfaDS[10] = -auxilliaryfXiM * auxilliaryfZetaP;
			auxilliaryfaDS[11] = -auxilliaryfXiP * auxilliaryfZetaP;
			auxilliaryfaDS[12] = -auxilliaryfaDS[15];
			auxilliaryfaDS[13] = -auxilliaryfaDS[14];
			auxilliaryfaDS[8] = -auxilliaryfaDS[11];
			auxilliaryfaDS[9] = -auxilliaryfaDS[10];


			auxilliaryfaDS[22] = -auxilliaryfXiM * auxilliaryfEtaM;
			auxilliaryfaDS[23] = -auxilliaryfXiP * auxilliaryfEtaM;
			auxilliaryfaDS[20] = -auxilliaryfXiP * auxilliaryfEtaP;
			auxilliaryfaDS[21] = -auxilliaryfXiM * auxilliaryfEtaP;
			auxilliaryfaDS[18] = -auxilliaryfaDS[22];
			auxilliaryfaDS[19] = -auxilliaryfaDS[23];
			auxilliaryfaDS[16] = -auxilliaryfaDS[20];
			auxilliaryfaDS[17] = -auxilliaryfaDS[21];

			return auxilliaryfaDS;
		}

		protected static double determinantTolerance = 0.00000001;
		//TODO: This should be handled by JacobianHexa8Reverse
		private Tuple<double[,], double[,], double> CalcH8JDetJ(double[,] faXYZ, double[] faDS)
		{
			double[,] auxilliaryfaJ = new double[3, 3];
			auxilliaryfaJ[0, 0] = faDS[0] * faXYZ[0, 0] + faDS[1] * faXYZ[0, 1] + faDS[2] * faXYZ[0, 2] + faDS[3] * faXYZ[0, 3] + faDS[4] * faXYZ[0, 4] + faDS[5] * faXYZ[0, 5] + faDS[6] * faXYZ[0, 6] + faDS[7] * faXYZ[0, 7];
			auxilliaryfaJ[0, 1] = faDS[0] * faXYZ[1, 0] + faDS[1] * faXYZ[1, 1] + faDS[2] * faXYZ[1, 2] + faDS[3] * faXYZ[1, 3] + faDS[4] * faXYZ[1, 4] + faDS[5] * faXYZ[1, 5] + faDS[6] * faXYZ[1, 6] + faDS[7] * faXYZ[1, 7];
			auxilliaryfaJ[0, 2] = faDS[0] * faXYZ[2, 0] + faDS[1] * faXYZ[2, 1] + faDS[2] * faXYZ[2, 2] + faDS[3] * faXYZ[2, 3] + faDS[4] * faXYZ[2, 4] + faDS[5] * faXYZ[2, 5] + faDS[6] * faXYZ[2, 6] + faDS[7] * faXYZ[2, 7];
			auxilliaryfaJ[1, 0] = faDS[8] * faXYZ[0, 0] + faDS[9] * faXYZ[0, 1] + faDS[10] * faXYZ[0, 2] + faDS[11] * faXYZ[0, 3] + faDS[12] * faXYZ[0, 4] + faDS[13] * faXYZ[0, 5] + faDS[14] * faXYZ[0, 6] + faDS[15] * faXYZ[0, 7];
			auxilliaryfaJ[1, 1] = faDS[8] * faXYZ[1, 0] + faDS[9] * faXYZ[1, 1] + faDS[10] * faXYZ[1, 2] + faDS[11] * faXYZ[1, 3] + faDS[12] * faXYZ[1, 4] + faDS[13] * faXYZ[1, 5] + faDS[14] * faXYZ[1, 6] + faDS[15] * faXYZ[1, 7];
			auxilliaryfaJ[1, 2] = faDS[8] * faXYZ[2, 0] + faDS[9] * faXYZ[2, 1] + faDS[10] * faXYZ[2, 2] + faDS[11] * faXYZ[2, 3] + faDS[12] * faXYZ[2, 4] + faDS[13] * faXYZ[2, 5] + faDS[14] * faXYZ[2, 6] + faDS[15] * faXYZ[2, 7];
			auxilliaryfaJ[2, 0] = faDS[16] * faXYZ[0, 0] + faDS[17] * faXYZ[0, 1] + faDS[18] * faXYZ[0, 2] + faDS[19] * faXYZ[0, 3] + faDS[20] * faXYZ[0, 4] + faDS[21] * faXYZ[0, 5] + faDS[22] * faXYZ[0, 6] + faDS[23] * faXYZ[0, 7];
			auxilliaryfaJ[2, 1] = faDS[16] * faXYZ[1, 0] + faDS[17] * faXYZ[1, 1] + faDS[18] * faXYZ[1, 2] + faDS[19] * faXYZ[1, 3] + faDS[20] * faXYZ[1, 4] + faDS[21] * faXYZ[1, 5] + faDS[22] * faXYZ[1, 6] + faDS[23] * faXYZ[1, 7];
			auxilliaryfaJ[2, 2] = faDS[16] * faXYZ[2, 0] + faDS[17] * faXYZ[2, 1] + faDS[18] * faXYZ[2, 2] + faDS[19] * faXYZ[2, 3] + faDS[20] * faXYZ[2, 4] + faDS[21] * faXYZ[2, 5] + faDS[22] * faXYZ[2, 6] + faDS[23] * faXYZ[2, 7];

			double auxilliaryfDet1 = auxilliaryfaJ[0, 0] * (auxilliaryfaJ[1, 1] * auxilliaryfaJ[2, 2] - auxilliaryfaJ[2, 1] * auxilliaryfaJ[1, 2]);
			double auxilliaryfDet2 = -auxilliaryfaJ[0, 1] * (auxilliaryfaJ[1, 0] * auxilliaryfaJ[2, 2] - auxilliaryfaJ[2, 0] * auxilliaryfaJ[1, 2]);
			double auxilliaryfDet3 = auxilliaryfaJ[0, 2] * (auxilliaryfaJ[1, 0] * auxilliaryfaJ[2, 1] - auxilliaryfaJ[2, 0] * auxilliaryfaJ[1, 1]);
			double auxilliaryfDetJ = auxilliaryfDet1 + auxilliaryfDet2 + auxilliaryfDet3;
			if (auxilliaryfDetJ < determinantTolerance)
			{
				throw new ArgumentException(
					$"Jacobian determinant is negative or under tolerance ({auxilliaryfDetJ} < {determinantTolerance})."
					 + " Check the order of nodes or the element geometry.");
			}

			double auxilliaryfDetInv = 1.0 / auxilliaryfDetJ;
			double[,] auxilliaryfaJInv = new double[3, 3];
			auxilliaryfaJInv[0, 0] = (auxilliaryfaJ[1, 1] * auxilliaryfaJ[2, 2] - auxilliaryfaJ[2, 1] * auxilliaryfaJ[1, 2]) * auxilliaryfDetInv;
			auxilliaryfaJInv[1, 0] = (auxilliaryfaJ[2, 0] * auxilliaryfaJ[1, 2] - auxilliaryfaJ[1, 0] * auxilliaryfaJ[2, 2]) * auxilliaryfDetInv;
			auxilliaryfaJInv[2, 0] = (auxilliaryfaJ[1, 0] * auxilliaryfaJ[2, 1] - auxilliaryfaJ[2, 0] * auxilliaryfaJ[1, 1]) * auxilliaryfDetInv;
			auxilliaryfaJInv[0, 1] = (auxilliaryfaJ[2, 1] * auxilliaryfaJ[0, 2] - auxilliaryfaJ[0, 1] * auxilliaryfaJ[2, 2]) * auxilliaryfDetInv;
			auxilliaryfaJInv[1, 1] = (auxilliaryfaJ[0, 0] * auxilliaryfaJ[2, 2] - auxilliaryfaJ[2, 0] * auxilliaryfaJ[0, 2]) * auxilliaryfDetInv;
			auxilliaryfaJInv[2, 1] = (auxilliaryfaJ[2, 0] * auxilliaryfaJ[0, 1] - auxilliaryfaJ[2, 1] * auxilliaryfaJ[0, 0]) * auxilliaryfDetInv;
			auxilliaryfaJInv[0, 2] = (auxilliaryfaJ[0, 1] * auxilliaryfaJ[1, 2] - auxilliaryfaJ[1, 1] * auxilliaryfaJ[0, 2]) * auxilliaryfDetInv;
			auxilliaryfaJInv[1, 2] = (auxilliaryfaJ[1, 0] * auxilliaryfaJ[0, 2] - auxilliaryfaJ[0, 0] * auxilliaryfaJ[1, 2]) * auxilliaryfDetInv;
			auxilliaryfaJInv[2, 2] = (auxilliaryfaJ[0, 0] * auxilliaryfaJ[1, 1] - auxilliaryfaJ[1, 0] * auxilliaryfaJ[0, 1]) * auxilliaryfDetInv;

			return new Tuple<double[,], double[,], double>(auxilliaryfaJ, auxilliaryfaJInv, auxilliaryfDetJ);
		}

		private readonly IList<INode> embeddedNodes = new List<INode>();
		public EmbeddedNode BuildHostElementEmbeddedNode(IElementType element, INode node,
			IEmbeddedDOFInHostTransformationVector transformationVector)
		{
			var points = GetNaturalCoordinates(element, node);
			if (points.Length == 0) return null;

			embeddedNodes.Add(node);
			var embeddedNode = new EmbeddedNode(node, element, transformationVector.GetDependentDOFTypes);
			for (int i = 0; i < points.Length; i++) embeddedNode.Coordinates.Add(points[i]);
			return embeddedNode;
		}

		private double[] GetNaturalCoordinates(IElementType element, INode node)
		{
			double[] mins = new double[] { element.Nodes[0].X, element.Nodes[0].Y, element.Nodes[0].Z };
			double[] maxes = new double[] { element.Nodes[0].X, element.Nodes[0].Y, element.Nodes[0].Z };
			for (int i = 0; i < element.Nodes.Count; i++)
			{
				mins[0] = mins[0] > element.Nodes[i].X ? element.Nodes[i].X : mins[0];
				mins[1] = mins[1] > element.Nodes[i].Y ? element.Nodes[i].Y : mins[1];
				mins[2] = mins[2] > element.Nodes[i].Z ? element.Nodes[i].Z : mins[2];
				maxes[0] = maxes[0] < element.Nodes[i].X ? element.Nodes[i].X : maxes[0];
				maxes[1] = maxes[1] < element.Nodes[i].Y ? element.Nodes[i].Y : maxes[1];
				maxes[2] = maxes[2] < element.Nodes[i].Z ? element.Nodes[i].Z : maxes[2];
			}
			//return new double[] { (node.X - mins[0]) / ((maxes[0] - mins[0]) / 2) - 1,
			//    (node.Y - mins[1]) / ((maxes[1] - mins[1]) / 2) - 1,
			//    (node.Z - mins[2]) / ((maxes[2] - mins[2]) / 2) - 1 };

			bool maybeInsideElement = node.X <= maxes[0] && node.X >= mins[0] &&
				node.Y <= maxes[1] && node.Y >= mins[1] &&
				node.Z <= maxes[2] && node.Z >= mins[2];
			if (maybeInsideElement == false) return new double[0];

			const int jacobianSize = 3;
			const int maxIterations = 1000;
			const double tolerance = 1e-10;
			int iterations = 0;
			double deltaNaturalCoordinatesNormSquare = 100;
			double[] naturalCoordinates = new double[] { 0, 0, 0 };
			const double toleranceSquare = tolerance * tolerance;

			while (deltaNaturalCoordinatesNormSquare > toleranceSquare && iterations < maxIterations)
			{
				iterations++;
				var shapeFunctions = CalcH8Shape(naturalCoordinates[0], naturalCoordinates[1], naturalCoordinates[2]);
				double[] coordinateDifferences = new double[] { 0, 0, 0 };
				for (int i = 0; i < shapeFunctions.Length; i++)
				{
					coordinateDifferences[0] += shapeFunctions[i] * element.Nodes[i].X;
					coordinateDifferences[1] += shapeFunctions[i] * element.Nodes[i].Y;
					coordinateDifferences[2] += shapeFunctions[i] * element.Nodes[i].Z;
				}
				coordinateDifferences[0] = node.X - coordinateDifferences[0];
				coordinateDifferences[1] = node.Y - coordinateDifferences[1];
				coordinateDifferences[2] = node.Z - coordinateDifferences[2];

				double[,] faXYZ = GetCoordinatesTranspose(element);
				double[] nablaShapeFunctions = CalcH8NablaShape(naturalCoordinates[0], naturalCoordinates[1], naturalCoordinates[2]);
				var inverseJacobian = CalcH8JDetJ(faXYZ, nablaShapeFunctions).Item2;

				double[] deltaNaturalCoordinates = new double[] { 0, 0, 0 };
				for (int i = 0; i < jacobianSize; i++)
					for (int j = 0; j < jacobianSize; j++)
						deltaNaturalCoordinates[i] += inverseJacobian[j, i] * coordinateDifferences[j];
				for (int i = 0; i < 3; i++)
					naturalCoordinates[i] += deltaNaturalCoordinates[i];

				deltaNaturalCoordinatesNormSquare = 0;
				for (int i = 0; i < 3; i++)
					deltaNaturalCoordinatesNormSquare += deltaNaturalCoordinates[i] * deltaNaturalCoordinates[i];
				//deltaNaturalCoordinatesNormSquare = Math.Sqrt(deltaNaturalCoordinatesNormSquare);
			}

			return naturalCoordinates.Count(x => Math.Abs(x) - 1.0 > tolerance) > 0 ? new double[0] : naturalCoordinates;
		}

		public double[] GetShapeFunctionsForNode(IElementType element, EmbeddedNode node)
		{
			double[,] elementCoordinates = GetCoordinatesTranspose(element);
			var shapeFunctions = CalcH8Shape(node.Coordinates[0], node.Coordinates[1], node.Coordinates[2]);
			var nablaShapeFunctions = CalcH8NablaShape(node.Coordinates[0], node.Coordinates[1], node.Coordinates[2]);
			var jacobian = CalcH8JDetJ(elementCoordinates, nablaShapeFunctions);

			return new double[]
			{
				shapeFunctions[0], shapeFunctions[1], shapeFunctions[2], shapeFunctions[3], shapeFunctions[4], shapeFunctions[5], shapeFunctions[6], shapeFunctions[7],
				nablaShapeFunctions[0], nablaShapeFunctions[1], nablaShapeFunctions[2], nablaShapeFunctions[3], nablaShapeFunctions[4], nablaShapeFunctions[5], nablaShapeFunctions[6], nablaShapeFunctions[7],
				nablaShapeFunctions[8], nablaShapeFunctions[9], nablaShapeFunctions[10], nablaShapeFunctions[11], nablaShapeFunctions[12], nablaShapeFunctions[13], nablaShapeFunctions[14], nablaShapeFunctions[15],
				nablaShapeFunctions[16], nablaShapeFunctions[17], nablaShapeFunctions[18], nablaShapeFunctions[19], nablaShapeFunctions[20], nablaShapeFunctions[21], nablaShapeFunctions[22], nablaShapeFunctions[23],
				jacobian.Item1[0, 0], jacobian.Item1[0, 1], jacobian.Item1[0, 2], jacobian.Item1[1, 0], jacobian.Item1[1, 1], jacobian.Item1[1, 2], jacobian.Item1[2, 0], jacobian.Item1[2, 1], jacobian.Item1[2, 2],
				jacobian.Item2[0, 0], jacobian.Item2[0, 1], jacobian.Item2[0, 2], jacobian.Item2[1, 0], jacobian.Item2[1, 1], jacobian.Item2[1, 2], jacobian.Item2[2, 0], jacobian.Item2[2, 1], jacobian.Item2[2, 2]
			};
		}

		#endregion


		public IEnumerable<IEnumerable<double>> InterpolateElementModelQuantities(IEnumerable<IElementModelQuantity<IStructuralDofType>> quantities, IEnumerable<double[]> coordinates) =>
			throw new NotImplementedException();
		public IEnumerable<double[]> IntegrateElementModelQuantities(IEnumerable<IElementModelQuantity<IStructuralDofType>> quantities) =>
			throw new NotImplementedException();
	}


}
