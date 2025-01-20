//TODO: Is there any point in having different material properties per Gauss point?
namespace MGroup.FEM.Thermal.Isoparametric
{
	using System;
	using System.Collections.Generic;

	using MGroup.Constitutive.Thermal;
	using MGroup.LinearAlgebra.Implementations;
	using MGroup.LinearAlgebra.Matrices;
	using MGroup.MSolve.DataStructures;
	using MGroup.MSolve.Discretization;
	using MGroup.MSolve.Discretization.Dofs;
	using MGroup.MSolve.Discretization.Embedding;
	using MGroup.MSolve.Discretization.Entities;
	using MGroup.MSolve.Geometry.Coordinates;
	using MGroup.MSolve.Numerics.Integration.Quadratures;
	using MGroup.MSolve.Numerics.Interpolation;
	using MGroup.MSolve.Numerics.Interpolation.GaussPointExtrapolation;
	using MGroup.MSolve.Numerics.Interpolation.Inverse;
	using MGroup.MSolve.Numerics.Interpolation.Jacobians;

	public class ThermalElement2D : IThermalElementType, IEmbeddedHostElement
	{
		private readonly IDofType[][] dofTypes; //TODO: this should not be stored for each element. Instead store it once for each Quad4, Tri3, etc. Otherwise create it on the fly.
		private readonly IThermalProperties material;
		//private readonly Dictionary<GaussPoint, ThermalMaterial> materialsAtGaussPoints;


		public ThermalElement2D(double thickness, IReadOnlyList<Node> nodes, IIsoparametricInterpolation2D interpolation,
			IQuadrature2D quadratureForStiffness, IQuadrature2D quadratureForConsistentMass,
			IGaussPointExtrapolation2D gaussPointExtrapolation,
			IThermalProperties material)
		{
			this.material = material;
			this.GaussPointExtrapolation = gaussPointExtrapolation;
			this.Nodes = nodes;
			this.Interpolation = interpolation;
			this.QuadratureForConsistentMass = quadratureForConsistentMass;
			this.QuadratureForStiffness = quadratureForStiffness;
			this.Thickness = thickness;

			dofTypes = new IDofType[nodes.Count][];
			for (int i = 0; i < interpolation.NumFunctions; ++i) dofTypes[i] = new IDofType[] { ThermalDof.Temperature };
		}

		public CellType CellType => Interpolation.CellType;
		public ElementDimensions ElementDimensions => ElementDimensions.TwoD;

		public int ID { get; set; }
		public int SubdomainID { get; set; }
		public IReadOnlyList<INode> Nodes { get; }

		public IGaussPointExtrapolation2D GaussPointExtrapolation { get; }
		public IIsoparametricInterpolation2D Interpolation { get; }
		public IQuadrature2D QuadratureForConsistentMass { get; }
		public IQuadrature2D QuadratureForStiffness { get; }
		public double Thickness { get; }

		public bool ConstitutiveLawModified => throw new NotImplementedException();

		public IElementDofEnumerator DofEnumerator { get; set; } = new GenericDofEnumerator();

		public IMatrix CapacityMatrix()
		{
			return BuildCapacityMatrix();
		}

		public Matrix BuildCapacityMatrix()
		{
			int numDofs = Nodes.Count;
			var capacity = Matrix.CreateZero(numDofs, numDofs);
			IReadOnlyList<double[]> shapeFunctions =
				Interpolation.EvaluateFunctionsAtGaussPoints(QuadratureForConsistentMass);
			IReadOnlyList<Matrix> shapeGradientsNatural =
				Interpolation.EvaluateNaturalGradientsAtGaussPoints(QuadratureForConsistentMass);

			for (int gp = 0; gp < QuadratureForConsistentMass.IntegrationPoints.Count; ++gp)
			{
				Matrix shapeFunctionMatrix = BuildShapeFunctionMatrix(shapeFunctions[gp]);
				Matrix partial = shapeFunctionMatrix.Transpose() * shapeFunctionMatrix;
				var jacobian = new IsoparametricJacobian2D(Nodes, shapeGradientsNatural[gp]);
				double dA = jacobian.DirectDeterminant * QuadratureForConsistentMass.IntegrationPoints[gp].Weight;
				capacity.AxpyIntoThis(partial, dA);
			}

			//WARNING: the following needs to change for non uniform density. Perhaps the integration order too.
			capacity.ScaleIntoThis(Thickness * material.Density * material.SpecialHeatCoeff);
			capacity.MatrixSymmetry = MatrixSymmetry.Symmetric;
			return capacity;
		}

		public Matrix BuildConductivityMatrix()
		{
			int numDofs = Nodes.Count;
			var conductivity = Matrix.CreateZero(numDofs, numDofs);
			IReadOnlyList<Matrix> shapeGradientsNatural =
				Interpolation.EvaluateNaturalGradientsAtGaussPoints(QuadratureForStiffness);

			for (int gp = 0; gp < QuadratureForStiffness.IntegrationPoints.Count; ++gp)
			{
				// Calculate the necessary quantities for the integration
				//Matrix constitutive = (Matrix)(materialsAtGaussPoints[gp].ConstitutiveMatrix); // ugly cast will be removed along with the legacy Matrix classes
				var jacobian = new IsoparametricJacobian2D(Nodes, shapeGradientsNatural[gp]);
				Matrix shapeGradientsCartesian =
					jacobian.TransformNaturalDerivativesToCartesian(shapeGradientsNatural[gp]);
				Matrix deformation = BuildDeformationMatrix(shapeGradientsCartesian);

				// Contribution of this gauss point to the element stiffness matrix
				Matrix partialK = deformation.Transpose() * deformation;
				//Matrix partialΚ = deformation.Transpose() * (constitutive * deformation);
				//partialK.Scale(materialsAtGaussPoints[gaussPoint].ThermalConductivity);

				double dA = jacobian.DirectDeterminant * QuadratureForStiffness.IntegrationPoints[gp].Weight; //TODO: this is used by all methods that integrate. I should cache it.
				conductivity.AxpyIntoThis(partialK, dA * material.ThermalConductivity);
			}

			conductivity.ScaleIntoThis(Thickness);
			conductivity.MatrixSymmetry = MatrixSymmetry.Symmetric;
			return conductivity;
		}

		// Provatidis uses two distinct vectors K = N,x^T * k * N,x + N,y^T * k * N,y
		//private (Matrix dNdX, Matrix dNdY) CalcdNdx(EvalShapeGradients2D shapeGrad)
		//{
		//    int n = Nodes.Count;
		//    var dNdX = new double[n, 1];
		//    var dNdY = new double[n, 1];
		//    for (int i = 0; i < n; ++i)
		//    {
		//        dNdX[i, 0] = shapeGrad[i][0];
		//        dNdY[i, 0] = shapeGrad[i][1];
		//    }
		//    return (new Matrix(dNdX), new Matrix(dNdY));
		//}

		private Matrix BuildDeformationMatrix(Matrix shapeGradientsCartesian)
		{
			//TODO: isn't this just the transpose of [dNi/dxj]?
			var deformation = Matrix.CreateZero(2, Nodes.Count);
			for (int nodeIdx = 0; nodeIdx < Nodes.Count; ++nodeIdx)
			{
				deformation[0, nodeIdx] = shapeGradientsCartesian[nodeIdx, 0];
				deformation[1, nodeIdx] = shapeGradientsCartesian[nodeIdx, 1];
			}
			return deformation;
		}

		/// <summary>
		/// The shape function matrix is 1-by-n, where n = is the number of shape functions.
		/// </summary>
		public Matrix BuildShapeFunctionMatrix(double[] shapeFunctions) //TODO: reconsider this. As it is, it just returns the shape functions in a Matrix
		{
			//var array2D = new double[1, shapeFunctions.Length];
			//for (int i = 0; i < shapeFunctions.Length; ++i)
			//{
			//    array2D[0, i] = shapeFunctions[i];
			//}
			//return Matrix.CreateFromArray(array2D);
			return Matrix.CreateFromArray(shapeFunctions, 1, shapeFunctions.Length);
		}

		public IReadOnlyList<IReadOnlyList<IDofType>> GetElementDofTypes() => dofTypes;

		public void ResetConstitutiveLawModified()
		{
			throw new NotImplementedException();
		}

		public Tuple<double[], double[]> CalculateResponse(double[] localDisplacements)
		{
			throw new NotImplementedException();
		}

		public double[] CalculateResponseIntegral()
		{
			throw new NotImplementedException();
		}

		public double[] CalculateResponseIntegralForLogging(double[] localDisplacements)
		{
			throw new NotImplementedException();
		}

		//public double[] CalculateAccelerationResponseIntegral(IElementType element, IList<MassAccelerationLoad> loads)
		//{
		//	throw new NotImplementedException();
		//}

		public void SaveConstitutiveLawState(IHaveState externalState)
		{
			throw new NotImplementedException();
		}

		public void ClearConstitutiveLawState()
		{
			throw new NotImplementedException();
		}

		public void ClearConstitutiveLawStresses()
		{
			throw new NotImplementedException();
		}

		public IMatrix ConductivityMatrix()
		{
			return BuildConductivityMatrix();
		}

		public IMatrix PhysicsMatrix()
		{
			return ConductivityMatrix();
		}
		private readonly IList<INode> embeddedNodes = new List<INode>();

		public EmbeddedNode BuildHostElementEmbeddedNode(IElementType element, INode node, IEmbeddedDOFInHostTransformationVector transformationVector)
		{
			IInverseInterpolation2D inverseInterpolation = Interpolation.CreateInverseMappingFor(Nodes);
			double[] naturalCoordinates = inverseInterpolation.TransformPointCartesianToNatural(new CartesianPoint(node.X, node.Y)).Coordinates;

			if (naturalCoordinates.Length == 0) return null;

			embeddedNodes.Add(node);
			var embeddedNode = new EmbeddedNode(node, element, transformationVector.GetDependentDOFTypes);
			for (int i = 0; i < naturalCoordinates.Length; i++)
				embeddedNode.Coordinates.Add(naturalCoordinates[i]);
			return embeddedNode;
		}

		public double[] GetShapeFunctionsForNode(IElementType element, EmbeddedNode node)
		{
			return Interpolation.EvaluateFunctionsAt(new NaturalPoint(node.Coordinates[0], node.Coordinates[1]));

			//TODO: This method originally returned an array containing the shape functions, shape function derivatives and entries of the inverse jacobian matrix.
			//      a) This is retarded and extremely difficult to work with. b) This array is used in Hexa8TranslationAndRotationTransformationVector. Cohesive embedding 
			//      and regular embedding in thermal elements only need the shape functions, so we will return these for now.

			//double[,] elementCoordinates = GetCoordinatesTranspose(element);
			//var shapeFunctions = CalcH8Shape(node.Coordinates[0], node.Coordinates[1], node.Coordinates[2]);
			//var nablaShapeFunctions = CalcH8NablaShape(node.Coordinates[0], node.Coordinates[1], node.Coordinates[2]);
			//var jacobian = CalcH8JDetJ(elementCoordinates, nablaShapeFunctions);

			//return new double[]
			//{
			//    shapeFunctions[0], shapeFunctions[1], shapeFunctions[2], shapeFunctions[3], shapeFunctions[4], shapeFunctions[5], shapeFunctions[6], shapeFunctions[7],
			//    nablaShapeFunctions[0], nablaShapeFunctions[1], nablaShapeFunctions[2], nablaShapeFunctions[3], nablaShapeFunctions[4], nablaShapeFunctions[5], nablaShapeFunctions[6], nablaShapeFunctions[7],
			//    nablaShapeFunctions[8], nablaShapeFunctions[9], nablaShapeFunctions[10], nablaShapeFunctions[11], nablaShapeFunctions[12], nablaShapeFunctions[13], nablaShapeFunctions[14], nablaShapeFunctions[15],
			//    nablaShapeFunctions[16], nablaShapeFunctions[17], nablaShapeFunctions[18], nablaShapeFunctions[19], nablaShapeFunctions[20], nablaShapeFunctions[21], nablaShapeFunctions[22], nablaShapeFunctions[23],
			//    jacobian.Item1[0, 0], jacobian.Item1[0, 1], jacobian.Item1[0, 2], jacobian.Item1[1, 0], jacobian.Item1[1, 1], jacobian.Item1[1, 2], jacobian.Item1[2, 0], jacobian.Item1[2, 1], jacobian.Item1[2, 2],
			//    jacobian.Item2[0, 0], jacobian.Item2[0, 1], jacobian.Item2[0, 2], jacobian.Item2[1, 0], jacobian.Item2[1, 1], jacobian.Item2[1, 2], jacobian.Item2[2, 0], jacobian.Item2[2, 1], jacobian.Item2[2, 2]
			//};
		}
	}
}
