using System;
using System.Collections.Generic;
using MGroup.Constitutive.Structural;
using MGroup.MSolve.Numerics.Interpolation;
using MGroup.MSolve.Numerics.Interpolation.Jacobians;
using MGroup.LinearAlgebra.Matrices;
using MGroup.Constitutive.Structural.Continuum;
using MGroup.MSolve.Discretization;
using MGroup.MSolve.Discretization.Dofs;
using MGroup.MSolve.Numerics.Integration.Quadratures;
using MGroup.MSolve.Discretization.BoundaryConditions;
using MGroup.LinearAlgebra.Vectors;
using System.Linq;
using MGroup.MSolve.Discretization.Entities;
using MGroup.FEM.Structural.Helpers;
using MGroup.MSolve.DataStructures;
using MGroup.MSolve.Constitutive;
using MGroup.LinearAlgebra.Implementations;

namespace MGroup.FEM.Structural.Continuum
{
	/// <summary>
	/// Continuum finite Element for 3d problems with material and geometric nonlinearities
	/// Authors: Gerasimos Sotiropoulos
	/// </summary>
	public class ContinuumElement3DGrowth : IStructuralElementType
	{
		protected readonly IDofType[] nodalDOFTypes = new IDofType[] { StructuralDof.TranslationX, StructuralDof.TranslationY, StructuralDof.TranslationZ };
		protected readonly IDofType[][] dofTypes;
		protected readonly IContinuumMaterial3DDefGrad[] materialsAtGaussPoints;
		protected IElementDofEnumerator dofEnumerator = new GenericDofEnumerator();

		private readonly int nGaussPoints;
		private bool isInitialized = false;

		private double[][] initialCoordinates; //not defined by user. 8 arrays of 3 elements
		private double[][] totalDisplacements;
		private double[] integrationCoeffs;

		private double[][] strainsVec;
		private double[][] lastStresses;
		private double[][] DefGradVec;
		private double lambdag = 1;

		public ContinuumElement3DGrowth(IReadOnlyList<INode> nodes, IContinuumMaterial3DDefGrad material, IQuadrature3D quadratureForStiffness,
			 IIsoparametricInterpolation3D interpolation, double lambdag = 1)
		{
			this.nGaussPoints = quadratureForStiffness.IntegrationPoints.Count;
			this.QuadratureForStiffness = quadratureForStiffness;
			this.Interpolation = interpolation;
			this.Nodes = nodes;
			this.lambdag = lambdag;

			lastStresses = new double[nGaussPoints][];
			materialsAtGaussPoints = new IContinuumMaterial3DDefGrad[nGaussPoints];
			for (int i = 0; i < nGaussPoints; i++)
			{
				materialsAtGaussPoints[i] = (IContinuumMaterial3DDefGrad)material.Clone();
				lastStresses[i] = new double[6];
			}

			dofTypes = new IDofType[nodes.Count][];
			for (int i = 0; i < nodes.Count; i++)
			{
				dofTypes[i] = new IDofType[]
				{
					StructuralDof.TranslationX, StructuralDof.TranslationY, StructuralDof.TranslationZ
				};
			}
		}

		private IIsoparametricInterpolation3D Interpolation { get; }
		private IQuadrature3D QuadratureForStiffness { get; }
		private IQuadrature3D QuadratureForConsistentMass { get; }
		public int ID { get; set; }
		public int SubdomainID { get; set; }
		public IReadOnlyList<INode> Nodes { get; }
		public CellType CellType => Interpolation.CellType;

		public ElementDimensions ElementDimensions => ElementDimensions.ThreeD;

		public IElementDofEnumerator DofEnumerator
		{
			get { return dofEnumerator; }
			set { dofEnumerator = value; }
		}

		//public IReadOnlyList<IStructuralMaterial> Materials => materialsAtGaussPoints;

		private Matrix[] Getbl13DeformationMatrices(IReadOnlyList<Matrix> shapeFunctionNaturalDerivatives)
		{
			Matrix[] bl13Matrices;
			bl13Matrices = new Matrix[nGaussPoints];
			for (int npoint = 0; npoint < nGaussPoints; npoint++)
			{
				bl13Matrices[npoint] = Matrix.CreateZero(9, 3 * shapeFunctionNaturalDerivatives[npoint].NumRows);
				for (int m = 0; m < shapeFunctionNaturalDerivatives[npoint].NumRows; m++)
				{
					for (int n = 0; n < 3; n++)
					{
						bl13Matrices[npoint][n, 3 * m + 0] = shapeFunctionNaturalDerivatives[npoint][m, n];
						bl13Matrices[npoint][n + 3, 3 * m + 1] = shapeFunctionNaturalDerivatives[npoint][m, n];
						bl13Matrices[npoint][n + 6, 3 * m + 2] = shapeFunctionNaturalDerivatives[npoint][m, n];
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
			int numNodes = Nodes.Count;
			IReadOnlyList<Matrix> shapeFunctionNaturalDerivatives;
			shapeFunctionNaturalDerivatives = Interpolation.EvaluateNaturalGradientsAtGaussPoints(QuadratureForStiffness);
			Matrix[] bl13Matrices;
			bl13Matrices = Getbl13DeformationMatrices(shapeFunctionNaturalDerivatives);

			Matrix[] bnl1Matrices;
			initialCoordinates = new double[numNodes][];
			totalDisplacements = new double[numNodes][];

			var jacobians = shapeFunctionNaturalDerivatives.Select(x => new IsoparametricJacobian3D(Nodes, x));

			Matrix[] jacobianInverse = jacobians.Select(x => x.InverseMatrix.Transpose()).ToArray();
			double[] jacobianDeterminants = jacobians.Select(x => x.DirectDeterminant).ToArray();

			integrationCoeffs = new double[nGaussPoints];
			bnl1Matrices = GetAuxilliaryDeformationbnl1Matrices(jacobianInverse);
			for (int j = 0; j < numNodes; j++)
			{
				initialCoordinates[j] = new double[] { Nodes[j].X, Nodes[j].Y, Nodes[j].Z, };
			}

			for (int gpoint = 0; gpoint < nGaussPoints; gpoint++)
			{
				integrationCoeffs[gpoint] = jacobianDeterminants[gpoint] * QuadratureForStiffness.IntegrationPoints[gpoint].Weight;

			}

			totalDisplacements = new double[numNodes][];
			DefGradVec = new double[nGaussPoints][];
			for (int gpoint = 0; gpoint < nGaussPoints; gpoint++)
			{
				//strainsVec[gpoint] = new double[6]; //MS
				//strainsVec_last_converged[gpoint] = new double[6];
				DefGradVec[gpoint] = new double[9];
			}

			for (int k = 0; k < numNodes; k++)
			{
				totalDisplacements[k] = new double[3];
			}

			isInitialized = true;
		}

		private void UpdateCoordinateData(double[] localdisplacements, out double[][] deformedCoordinates)
		{
			int numNodes = localdisplacements.Length / 3;
			deformedCoordinates = new double[numNodes][];
			for (int j = 0; j < numNodes; j++)
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
			var jacobians = shapeFunctionNaturalDerivatives.Select(x => new IsoparametricJacobian3D(Nodes, x));
			Matrix[] jacobianInverse = jacobians.Select(x => x.InverseMatrix.Transpose()).ToArray();
			double[] jacobianDeterminants = jacobians.Select(x => x.DirectDeterminant).ToArray();
			//TODO: possibility of caching shapeFunctionNaturalDerivatives or J_0inv

			Matrix[] deformationGradientsTransposed = new Matrix[nGaussPoints];
			//Matrix[] GL = new Matrix[nGaussPoints];
			for (int npoint = 0; npoint < nGaussPoints; npoint++)
			{
				deformationGradientsTransposed[npoint] = Matrix.CreateZero(3, 3);
			}

			var jacobiansDeformed = shapeFunctionNaturalDerivatives.Select(x => new IsoparametricJacobian3D(deformedCoordinates, x, false)).ToArray();
			Matrix[] jacobiansDeformedMatrices = jacobiansDeformed.Select(x => x.DirectMatrix).ToArray();

			for (int npoint = 0; npoint < nGaussPoints; npoint++)
			{
				deformationGradientsTransposed[npoint] = jacobianInverse[npoint] * jacobiansDeformedMatrices[npoint];
				DefGradVec[npoint] = new double[9] { deformationGradientsTransposed[npoint][0, 0], deformationGradientsTransposed[npoint][1, 1],
					deformationGradientsTransposed[npoint][2, 2], deformationGradientsTransposed[npoint][1, 0], deformationGradientsTransposed[npoint][2, 1],
					deformationGradientsTransposed[npoint][0, 2], deformationGradientsTransposed[npoint][2, 0], deformationGradientsTransposed[npoint][0, 1],
					deformationGradientsTransposed[npoint][1, 2], };//MS
			}
		}

		private double[] UpdateResponseIntegral()
		{
			//TODO: the gauss point loop should be the outer one
			// Matrices that are not currently cached are calculated here.
			int numNodes = Nodes.Count();
			Matrix ll2 = Matrix.CreateZero(numNodes, 3);
			for (int m = 0; m < numNodes; m++)
			{
				for (int n = 0; n < 3; n++)
				{
					ll2[m, n] = totalDisplacements[m][n];
				}
			}

			IReadOnlyList<Matrix> shapeFunctionNaturalDerivatives;
			shapeFunctionNaturalDerivatives = Interpolation.EvaluateNaturalGradientsAtGaussPoints(QuadratureForStiffness);
			var jacobians = shapeFunctionNaturalDerivatives.Select(x => new IsoparametricJacobian3D(Nodes, x));
			Matrix[] jacobianInverse = jacobians.Select(x => x.InverseMatrix.Transpose()).ToArray();
			double[] jacobianDeterminants = jacobians.Select(x => x.DirectDeterminant).ToArray();

			Matrix[] bl13Matrices = Getbl13DeformationMatrices(shapeFunctionNaturalDerivatives);
			Matrix[] bl11aMatrices; // dimension number of gpoints
			Matrix[] bl12Marices;
			Matrix[] bl01Matrices;
			bl11aMatrices = Getbl11aDeformationMatrices(jacobianInverse);
			bl12Marices = GetBL12DeformationMatrices(jacobianInverse);
			bl01Matrices = Getbl01MDeformationMatrices(jacobianInverse);

			//INITIALIZATION of matrices that are currently not cached
			double[][] integrCoeffsTimesStresses = new double[nGaussPoints][];
			Matrix[] blMatrices = new Matrix[nGaussPoints];
			for (int gpoint = 0; gpoint < nGaussPoints; gpoint++)
			{
				integrCoeffsTimesStresses[gpoint] = new double[6];
				blMatrices[gpoint] = Matrix.CreateZero(6, 3 * numNodes);
			}

			double[][] forces = new double[nGaussPoints + 1][];
			for (int npoint = 0; npoint < nGaussPoints + 1; npoint++)
			{
				forces[npoint] = new double[3 * numNodes];
			}

			Matrix[] bl11Matrices = new Matrix[nGaussPoints];
			Matrix[] bL1112Plus01Matrices = new Matrix[nGaussPoints];
			for (int npoint = 0; npoint < nGaussPoints; npoint++)
			{
				bl11Matrices[npoint] = Matrix.CreateZero(6, 9);
				bL1112Plus01Matrices[npoint] = Matrix.CreateZero(6, 9);
			}

			for (int npoint = 0; npoint < nGaussPoints; npoint++)
			{
				var stressesElastic = lastStresses[npoint];
				var secondPiolaElastic = new double[3, 3] { {stressesElastic[0],stressesElastic[3],stressesElastic[5] },
					{ stressesElastic[3],stressesElastic[1],stressesElastic[4] },
					{stressesElastic[5],stressesElastic[4],stressesElastic[2] } };
				var defGradTransposed = new double[3, 3] {{DefGradVec[npoint][0],
				DefGradVec[npoint][7], DefGradVec[npoint][5] },
					{ DefGradVec[npoint][3], DefGradVec[npoint][1],
					DefGradVec[npoint][8] }, {DefGradVec[npoint][6],
					DefGradVec[npoint][4], DefGradVec[npoint][2] }, };
				var defGradElasticTransposed = new double[3, 3] {{DefGradVec[npoint][0] / lambdag,
				DefGradVec[npoint][7] / lambdag, DefGradVec[npoint][5] / lambdag },
					{ DefGradVec[npoint][3] / lambdag, DefGradVec[npoint][1] / lambdag,
					DefGradVec[npoint][8] / lambdag}, {DefGradVec[npoint][6] / lambdag,
					DefGradVec[npoint][4] / lambdag, DefGradVec[npoint][2] / lambdag }, };
				var firstPiolaElastic = new double[3, 3];
				for (int i = 0; i < 3; i++)
				{
					for (int j = 0; j < 3; j++)
					{
						for (int k = 0; k < 3; k++)
						{
							firstPiolaElastic[i, j] += secondPiolaElastic[i, k] * defGradElasticTransposed[k, j];
						}
					}
				}
				double defGradDeterminant = 0;
				for (int i = 0; i < 3; i++)
					defGradDeterminant = defGradDeterminant + (defGradTransposed[0, i] * (defGradTransposed[1, (i + 1) % 3]
						* defGradTransposed[2, (i + 2) % 3] - defGradTransposed[1, (i + 2) % 3] * defGradTransposed[2, (i + 1) % 3]));

				double defGradElasticDeterminant = 0;
				for (int i = 0; i < 3; i++)
					defGradElasticDeterminant = defGradElasticDeterminant + (defGradElasticTransposed[0, i] * (defGradElasticTransposed[1, (i + 1) % 3]
						* defGradElasticTransposed[2, (i + 2) % 3] - defGradElasticTransposed[1, (i + 2) % 3] * defGradElasticTransposed[2, (i + 1) % 3]));

				var detFg = Math.Pow(lambdag, 3);
				double[,] defGrad_lamdag_inv = new double[3, 3] { { (double)1 / lambdag, 0, 0 }, { 0, (double)1 / lambdag, 0 }, { 0, 0, (double)1 / lambdag } };

				var firstPiola = new double[3, 3];
				for (int i = 0; i < 3; i++)
				{
					for (int j = 0; j < 3; j++)
					{
						for (int k = 0; k < 3; k++)
						{
							firstPiola[i, j] += firstPiolaElastic[i, k] * defGrad_lamdag_inv[k, j] * detFg;
						}
					}
				}

				double[,] defGradInverseTransposed = new double[3, 3] { { (defGradTransposed[2, 2]*defGradTransposed[1,1]- defGradTransposed[2, 1]
					* defGradTransposed[1, 2])/defGradDeterminant,
						(-(defGradTransposed[2, 2] * defGradTransposed[0, 1] - defGradTransposed[2, 1]
					* defGradTransposed[0, 2]))/defGradDeterminant,
						(defGradTransposed[1,2] * defGradTransposed[0, 1] - defGradTransposed[1, 1] * defGradTransposed[0, 2])/defGradDeterminant },
					{ (-(defGradTransposed[2,2]*defGradTransposed[1,0]-defGradTransposed[2,0]*defGradTransposed[1,2]))/defGradDeterminant,
						(defGradTransposed[2,2]*defGradTransposed[0,0]-defGradTransposed[2,0]*defGradTransposed[0,2])/ defGradDeterminant,
						(-(defGradTransposed[1,2]*defGradTransposed[0,0]-defGradTransposed[1,0]*defGradTransposed[0,2]))/defGradDeterminant },
					{(defGradTransposed[2,1]*defGradTransposed[1,0]-defGradTransposed[2,0]*defGradTransposed[1,1])/defGradDeterminant,
						(-(defGradTransposed[2,1]*defGradTransposed[0,0]-defGradTransposed[2,0]*defGradTransposed[0,1]))/defGradDeterminant,
						(defGradTransposed[1,1]*defGradTransposed[0,0]-defGradTransposed[1,0]*defGradTransposed[0,1])/ defGradDeterminant } };

				var secondPiola = new double[3, 3];
				for (int i = 0; i < 3; i++)
				{
					for (int j = 0; j < 3; j++)
					{
						for (int k = 0; k < 3; k++)
						{
							secondPiola[i, j] += firstPiola[i, k] * defGradInverseTransposed[k, j];
						}
					}
				}

				var stresses = new double[6] { secondPiola[0, 0], secondPiola[1, 1], secondPiola[2, 2], secondPiola[0, 1], secondPiola[1, 2], secondPiola[2, 0] };
				integrCoeffsTimesStresses[npoint] = stresses.Scale(integrationCoeffs[npoint]);

				Matrix lcyrcumflex;//= Matrix.CreateZero(3, 3);
				lcyrcumflex = shapeFunctionNaturalDerivatives[npoint].Transpose() * ll2;

				for (int m = 0; m < 6; m++)
				{
					for (int n = 0; n < 3; n++)
					{
						for (int p = 0; p < 3; p++)
						{
							bl11Matrices[npoint][m, n] += bl11aMatrices[npoint][m, p] * lcyrcumflex[p, n];
							bl11Matrices[npoint][m, 3 + n] += bl11aMatrices[npoint][m, 3 + p] * lcyrcumflex[p, n];
							bl11Matrices[npoint][m, 6 + n] += bl11aMatrices[npoint][m, 6 + p] * lcyrcumflex[p, n];
						}
					}
				}

				bL1112Plus01Matrices[npoint] = bl11Matrices[npoint] * bl12Marices[npoint];
				bL1112Plus01Matrices[npoint].AddIntoThis(bl01Matrices[npoint]);

				blMatrices[npoint] = bL1112Plus01Matrices[npoint] * bl13Matrices[npoint];

				forces[npoint] = blMatrices[npoint].Multiply(integrCoeffsTimesStresses[npoint], true);
			}

			for (int npoint = 0; npoint < nGaussPoints; npoint++)
			{
				forces[nGaussPoints].AddIntoThis(forces[npoint]);
			}

			return forces[nGaussPoints];
		}

		private Matrix UpdateKmatrices()
		{
			int numNodes = Nodes.Count();
			Matrix elementStiffnessMatrix = Matrix.CreateZero(3 * numNodes, 3 * numNodes);


			// initialization of matrices that are not cached currently
			double[][] integrCoeffsTimesSpkvec = new double[nGaussPoints][];
			Matrix[] blMatrices = new Matrix[nGaussPoints];
			for (int gpoint = 0; gpoint < nGaussPoints; gpoint++)
			{
				integrCoeffsTimesSpkvec[gpoint] = new double[6];
				blMatrices[gpoint] = Matrix.CreateZero(6, 3 * numNodes);

			}
			Matrix totalDisplacementsMatrixReordered = Matrix.CreateZero(numNodes, 3);
			for (int m = 0; m < numNodes; m++)
			{
				for (int n = 0; n < 3; n++)
				{
					totalDisplacementsMatrixReordered[m, n] = totalDisplacements[m][n];
				}
			}
			IReadOnlyList<Matrix> shapeFunctionNaturalDerivatives;
			shapeFunctionNaturalDerivatives = Interpolation.EvaluateNaturalGradientsAtGaussPoints(QuadratureForStiffness);
			var jacobians = shapeFunctionNaturalDerivatives.Select(x => new IsoparametricJacobian3D(Nodes, x));
			Matrix[] jacobianInverse = jacobians.Select(x => x.InverseMatrix.Transpose()).ToArray();
			double[] jacobianDeterminants = jacobians.Select(x => x.DirectDeterminant).ToArray();

			Matrix[] bl13Matrices;
			bl13Matrices = Getbl13DeformationMatrices(shapeFunctionNaturalDerivatives);
			Matrix[] bl11aMatrices; // dimension: gpoints
			Matrix[] bl12Marices;
			Matrix[] bl01Matrices;
			bl11aMatrices = Getbl11aDeformationMatrices(jacobianInverse);
			bl12Marices = GetBL12DeformationMatrices(jacobianInverse);
			bl01Matrices = Getbl01MDeformationMatrices(jacobianInverse);

			Matrix[] bl11Matrices = new Matrix[nGaussPoints];
			Matrix[] bL1112Plus01Matrices = new Matrix[nGaussPoints];
			for (int npoint = 0; npoint < nGaussPoints; npoint++)
			{
				bl11Matrices[npoint] = Matrix.CreateZero(6, 9);
				bL1112Plus01Matrices[npoint] = Matrix.CreateZero(6, 9); //TODO this may be unnescessary
			}


			var secPiolaMat = new double[nGaussPoints][,];
			var secPiolaElasMat = new double[nGaussPoints][,];
			var s = MatrixSymmetry.Symmetric;

			for (int npoint = 0; npoint < nGaussPoints; npoint++)
			{
				var stressesElastic = lastStresses[npoint];
				secPiolaElasMat[npoint] = new double[3, 3] { {stressesElastic[0],stressesElastic[3],stressesElastic[5] },
					{ stressesElastic[3],stressesElastic[1],stressesElastic[4] },
					{stressesElastic[5],stressesElastic[4],stressesElastic[2] } };
				var defGradTransposed = new double[3, 3] {{DefGradVec[npoint][0],
				DefGradVec[npoint][7], DefGradVec[npoint][5] },
					{ DefGradVec[npoint][3], DefGradVec[npoint][1],
					DefGradVec[npoint][8] }, {DefGradVec[npoint][6],
					DefGradVec[npoint][4], DefGradVec[npoint][2] }, };
				var defGradElasticTransposed = new double[3, 3] {{DefGradVec[npoint][0] / lambdag,
				DefGradVec[npoint][7] / lambdag, DefGradVec[npoint][5] / lambdag },
					{ DefGradVec[npoint][3] / lambdag, DefGradVec[npoint][1] / lambdag,
					DefGradVec[npoint][8] / lambdag}, {DefGradVec[npoint][6] / lambdag,
					DefGradVec[npoint][4] / lambdag, DefGradVec[npoint][2] / lambdag }, };
				var firstPiolaElastic = new double[3, 3];
				for (int i = 0; i < 3; i++)
				{
					for (int j = 0; j < 3; j++)
					{
						for (int k = 0; k < 3; k++)
						{
							firstPiolaElastic[i, j] += secPiolaElasMat[npoint][i, k] * defGradElasticTransposed[k, j];
						}
					}
				}

				double defGradDeterminant = 0;
				for (int i = 0; i < 3; i++)
					defGradDeterminant = defGradDeterminant + (defGradTransposed[0, i] * (defGradTransposed[1, (i + 1) % 3]
						* defGradTransposed[2, (i + 2) % 3] - defGradTransposed[1, (i + 2) % 3] * defGradTransposed[2, (i + 1) % 3]));

				double defGradElasticDeterminant = 0;
				for (int i = 0; i < 3; i++)
					defGradElasticDeterminant = defGradElasticDeterminant + (defGradElasticTransposed[0, i] * (defGradElasticTransposed[1, (i + 1) % 3]
						* defGradElasticTransposed[2, (i + 2) % 3] - defGradElasticTransposed[1, (i + 2) % 3] * defGradElasticTransposed[2, (i + 1) % 3]));

				var detFg = Math.Pow(lambdag, 3);

				double[,] defGrad_lamdag_inv = new double[3, 3] { { (double)1 / lambdag, 0, 0 }, { 0, (double)1 / lambdag, 0 }, { 0, 0, (double)1 / lambdag } };

				var firstPiola = new double[3, 3];
				for (int i = 0; i < 3; i++)
				{
					for (int j = 0; j < 3; j++)
					{
						for (int k = 0; k < 3; k++)
						{
							firstPiola[i, j] += firstPiolaElastic[i, k] * defGrad_lamdag_inv[k, j] * detFg;
						}
					}
				}

				double[,] defGradInverseTransposed = new double[3, 3] { { (defGradTransposed[2, 2]*defGradTransposed[1,1]- defGradTransposed[2, 1]
					* defGradTransposed[1, 2])/defGradDeterminant,
						(-(defGradTransposed[2, 2] * defGradTransposed[0, 1] - defGradTransposed[2, 1]
					* defGradTransposed[0, 2]))/defGradDeterminant,
						(defGradTransposed[1,2] * defGradTransposed[0, 1] - defGradTransposed[1, 1] * defGradTransposed[0, 2])/defGradDeterminant },
					{ (-(defGradTransposed[2,2]*defGradTransposed[1,0]-defGradTransposed[2,0]*defGradTransposed[1,2]))/defGradDeterminant,
						(defGradTransposed[2,2]*defGradTransposed[0,0]-defGradTransposed[2,0]*defGradTransposed[0,2])/ defGradDeterminant,
						(-(defGradTransposed[1,2]*defGradTransposed[0,0]-defGradTransposed[1,0]*defGradTransposed[0,2]))/defGradDeterminant },
					{(defGradTransposed[2,1]*defGradTransposed[1,0]-defGradTransposed[2,0]*defGradTransposed[1,1])/defGradDeterminant,
						(-(defGradTransposed[2,1]*defGradTransposed[0,0]-defGradTransposed[2,0]*defGradTransposed[0,1]))/defGradDeterminant,
						(defGradTransposed[1,1]*defGradTransposed[0,0]-defGradTransposed[1,0]*defGradTransposed[0,1])/ defGradDeterminant } };

				secPiolaMat[npoint] = new double[3, 3];
				for (int i = 0; i < 3; i++)
				{
					for (int j = 0; j < 3; j++)
					{
						for (int k = 0; k < 3; k++)
						{
							secPiolaMat[npoint][i, j] += firstPiola[i, k] * defGradInverseTransposed[k, j];
						}
					}
				}

				var stresses = new double[6] { secPiolaMat[npoint][0, 0], secPiolaMat[npoint][1, 1], secPiolaMat[npoint][2, 2], secPiolaMat[npoint][0, 1], secPiolaMat[npoint][1, 2], secPiolaMat[npoint][2, 0] };

				integrCoeffsTimesSpkvec[npoint] = stresses.Scale(integrationCoeffs[npoint]);

				Matrix lcyrcumflex = Matrix.CreateZero(3, 3);
				lcyrcumflex = shapeFunctionNaturalDerivatives[npoint].Transpose() * totalDisplacementsMatrixReordered;

				for (int m = 0; m < 6; m++)
				{
					for (int n = 0; n < 3; n++)
					{
						for (int p = 0; p < 3; p++)
						{
							bl11Matrices[npoint][m, n] += bl11aMatrices[npoint][m, p] * lcyrcumflex[p, n];
							bl11Matrices[npoint][m, 3 + n] += bl11aMatrices[npoint][m, 3 + p] * lcyrcumflex[p, n];
							bl11Matrices[npoint][m, 6 + n] += bl11aMatrices[npoint][m, 6 + p] * lcyrcumflex[p, n];
						}
					}
				}

				bL1112Plus01Matrices[npoint] = bl11Matrices[npoint] * bl12Marices[npoint];
				bL1112Plus01Matrices[npoint].AddIntoThis(bl01Matrices[npoint]);

				blMatrices[npoint] = bL1112Plus01Matrices[npoint] * bl13Matrices[npoint];

			}
			// TODO: BL and above calculations can cached from calculate forces method

			Matrix[] bnl1Matrices;
			Matrix[] bnlMatrices;
			bnl1Matrices = GetAuxilliaryDeformationbnl1Matrices(jacobianInverse);
			bnlMatrices = new Matrix[nGaussPoints];
			for (int gpoint = 0; gpoint < nGaussPoints; gpoint++)
			{
				//bnlMatrices[gpoint] = Matrix.CreateZero(9, 3*numNodes); //todo this may be unnescessary

				bnlMatrices[gpoint] = bnl1Matrices[gpoint] * bl13Matrices[gpoint];

			}

			Matrix[] integrCoeffsTimesStresses = new Matrix[nGaussPoints];
			for (int npoint = 0; npoint < nGaussPoints; npoint++)
			{
				integrCoeffsTimesStresses[npoint] = Matrix.CreateZero(3, 3);
			}

			Matrix[] klStiffnessMatrixContributions = new Matrix[nGaussPoints + 1];
			Matrix[] knlStiffnessMatrixContributions = new Matrix[nGaussPoints + 1];
			for (int npoint = 0; npoint < nGaussPoints + 1; npoint++)
			{
				klStiffnessMatrixContributions[npoint] = Matrix.CreateZero(3 * numNodes, 3 * numNodes);
				knlStiffnessMatrixContributions[npoint] = Matrix.CreateZero(3 * numNodes, 3 * numNodes);
			}

			for (int npoint = 0; npoint < nGaussPoints; npoint++)
			{


				Matrix integrCoeffsTimesStressesTimesbnlMatrices = Matrix.CreateZero(9, 3 * numNodes); //TODO
				Matrix integrCoeffsTimesConsMatrix = Matrix.CreateZero(6, 6); //TODO
				Matrix integrCoeffTimesConsMatrixTimesBLMatrices = Matrix.CreateZero(6, 3 * numNodes);//TODO

				//
				integrCoeffsTimesStresses[npoint][0, 0] = integrCoeffsTimesSpkvec[npoint][0];
				integrCoeffsTimesStresses[npoint][0, 1] = integrCoeffsTimesSpkvec[npoint][3];
				integrCoeffsTimesStresses[npoint][0, 2] = integrCoeffsTimesSpkvec[npoint][5];
				integrCoeffsTimesStresses[npoint][1, 0] = integrCoeffsTimesSpkvec[npoint][3];
				integrCoeffsTimesStresses[npoint][1, 1] = integrCoeffsTimesSpkvec[npoint][1];
				integrCoeffsTimesStresses[npoint][1, 2] = integrCoeffsTimesSpkvec[npoint][4];
				integrCoeffsTimesStresses[npoint][2, 0] = integrCoeffsTimesSpkvec[npoint][5];
				integrCoeffsTimesStresses[npoint][2, 1] = integrCoeffsTimesSpkvec[npoint][4];
				integrCoeffsTimesStresses[npoint][2, 2] = integrCoeffsTimesSpkvec[npoint][2];

				//
				IMatrixView consDisp = materialsAtGaussPoints[npoint].ConstitutiveMatrix;
				s = s == MatrixSymmetry.Symmetric ? consDisp.MatrixSymmetry : s;

				var DG = new double[3, 3] { { DefGradVec[npoint][0], DefGradVec[npoint][3], DefGradVec[npoint][6] }, { DefGradVec[npoint][7], DefGradVec[npoint][1], DefGradVec[npoint][4] }, { DefGradVec[npoint][5], DefGradVec[npoint][8], DefGradVec[npoint][2] } };
				var DG_el = new double[3, 3] { { DefGradVec[npoint][0] / lambdag, DefGradVec[npoint][3] / lambdag, DefGradVec[npoint][6] / lambdag }, { DefGradVec[npoint][7] / lambdag, DefGradVec[npoint][1] / lambdag, DefGradVec[npoint][4] / lambdag }, { DefGradVec[npoint][5] / lambdag, DefGradVec[npoint][8] / lambdag, DefGradVec[npoint][2] / lambdag } };
				var Fg_inv = new double[3, 3] { { (double)1 / lambdag, 0, 0 }, { 0, (double)1 / lambdag, 0 }, { 0, 0, (double)1 / lambdag } };
				var detFg = Math.Pow(lambdag, 3);
				var consDisp2 = TransformationMethods.Calculate_dSdE_from_dSdEe(consDisp.CopytoArray2D(), DG, lambdag, secPiolaMat[npoint], DG_el, secPiolaElasMat[npoint], Fg_inv, detFg);
				consDisp = Matrix.CreateFromArray(consDisp2);

				for (int m = 0; m < 6; m++)
				{
					for (int n = 0; n < 6; n++)
					{
						integrCoeffsTimesConsMatrix[m, n] = integrationCoeffs[npoint] * consDisp[m, n];
					}
				}
				integrCoeffTimesConsMatrixTimesBLMatrices = integrCoeffsTimesConsMatrix * blMatrices[npoint];

				klStiffnessMatrixContributions[npoint] = blMatrices[npoint].Transpose() * integrCoeffTimesConsMatrixTimesBLMatrices;

				for (int m = 0; m < 3; m++) // 3x24 dimensions
				{
					for (int n = 0; n < 3 * numNodes; n++)
					{
						for (int p = 0; p < 3; p++)
						{
							integrCoeffsTimesStressesTimesbnlMatrices[m, n] += integrCoeffsTimesStresses[npoint][m, p] * bnlMatrices[npoint][p, n];
							integrCoeffsTimesStressesTimesbnlMatrices[3 + m, n] += integrCoeffsTimesStresses[npoint][m, p] * bnlMatrices[npoint][3 + p, n];
							integrCoeffsTimesStressesTimesbnlMatrices[6 + m, n] += integrCoeffsTimesStresses[npoint][m, p] * bnlMatrices[npoint][6 + p, n];
						}
					}
				}

				knlStiffnessMatrixContributions[npoint] = bnlMatrices[npoint].Transpose() * integrCoeffsTimesStressesTimesbnlMatrices;
			}

			// Add contributions of each gp on the total element stiffness matrix elementStiffnessMatrix            
			for (int npoint = 0; npoint < nGaussPoints; npoint++)
			{
				for (int m = 0; m < 3 * numNodes; m++)
				{
					for (int n = 0; n < 3 * numNodes; n++)
					{
						klStiffnessMatrixContributions[nGaussPoints][m, n] += klStiffnessMatrixContributions[npoint][m, n];
						knlStiffnessMatrixContributions[nGaussPoints][m, n] += knlStiffnessMatrixContributions[npoint][m, n];
					}
				}
			}

			for (int m = 0; m < 3 * numNodes; m++)
			{
				for (int n = 0; n < 3 * numNodes; n++)
				{
					elementStiffnessMatrix[m, n] = klStiffnessMatrixContributions[nGaussPoints][m, n] + knlStiffnessMatrixContributions[nGaussPoints][m, n];
				}
			}

			elementStiffnessMatrix.MatrixSymmetry = s;
			return elementStiffnessMatrix;
		}

		public Tuple<double[], double[]> CalculateResponse(double[] localTotalDisplacements)
		{
			this.UpdateCoordinateData(localTotalDisplacements, out double[][] deformedCoordinates);
			this.CalculateStrains(localTotalDisplacements, deformedCoordinates);
			for (int npoint = 0; npoint < materialsAtGaussPoints.Length; npoint++)
			{
				double[] DefGradVecEl = new double[9];
				for (int i = 0; i < 9; i++)
				{
					DefGradVecEl[i] = DefGradVec[npoint][i] / lambdag;
				}
				lastStresses[npoint] = materialsAtGaussPoints[npoint].UpdateConstitutiveMatrixAndEvaluateResponse(DefGradVecEl); //MS
			};

			return new Tuple<double[], double[]>(DefGradVec[materialsAtGaussPoints.Length - 1], lastStresses[materialsAtGaussPoints.Length - 1]);

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
				int numNodes = Nodes.Count();
				this.CalculateInitialConfigurationData();
				var localTotalDisplacements = new double[3 * numNodes];
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

		public void ClearConstitutiveLawState()
		{
			//TODO: the next throws an exception. Investigate. Possible changes in Analyzers may be the cause.
			//foreach (IContinuumMaterial3DDefGrad m in materialsAtGaussPoints) m.ClearState();
		}

		

		public void SaveConstitutiveLawState(IHaveState externalState)
		{
			//for (int npoint = 0; npoint < materialsAtGaussPoints.Length; npoint++)
			//{
			//	for (int i1 = 0; i1 < 6; i1++)
			//	{ strainsVecLastConverged[npoint][i1] = strainsVec[npoint][i1]; }
			//}

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
		//    foreach (IContinuumMaterial3D m in materialsAtGaussPoints) m.ClearStresses();
		//}

		public virtual IReadOnlyList<IReadOnlyList<IDofType>> GetElementDofTypes() => dofTypes;

		#region not implemented
		//public double[] CalculateAccelerationResponseIntegral(IElementType element, IList<MassAccelerationLoad> loads)
		//{
		//    throw new NotImplementedException();
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

		public IEnumerable<IEnumerable<double>> InterpolateElementModelQuantities(IEnumerable<IElementModelQuantity<IStructuralDofType>> quantities, IEnumerable<double[]> coordinates) =>
			throw new NotImplementedException();
		public IEnumerable<double[]> IntegrateElementModelQuantities(IEnumerable<IElementModelQuantity<IStructuralDofType>> quantities) =>
			throw new NotImplementedException();
	}
}
