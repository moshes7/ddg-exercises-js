"use strict";

/**
 * This class contains methods to build common {@link https://cs.cmu.edu/~kmcrane/Projects/DDG/paper.pdf discrete exterior calculus} operators.
 * @memberof module:Core
 */
class DEC {
	/**
	 * Builds a sparse diagonal matrix encoding the Hodge operator on 0-forms.
	 * By convention, the area of a vertex is 1.
	 * @static
	 * @param {module:Core.Geometry} geometry The geometry of a mesh.
	 * @param {Object} vertexIndex A dictionary mapping each vertex of a mesh to a unique index.
	 * @returns {module:LinearAlgebra.SparseMatrix}
	 */
	static buildHodgeStar0Form(geometry, vertexIndex) {

		let V = geometry.mesh.vertices.length;
		let T = new Triplet(V, V);
		let i_v, v, A_dual;
		
		for (let i=0; i<V; i++) {
			i_v = vertexIndex[i];
			v = geometry.mesh.vertices[i_v];
			A_dual = geometry.barycentricDualArea(v);
			T.addEntry(A_dual, i_v, i_v); 
		}

		let hodge0Form = SparseMatrix.fromTriplet(T);
		
		return hodge0Form;
	}

	/**
	 * Builds a sparse diagonal matrix encoding the Hodge operator on 1-forms.
	 * @static
	 * @param {module:Core.Geometry} geometry The geometry of a mesh.
	 * @param {Object} edgeIndex A dictionary mapping each edge of a mesh to a unique index.
	 * @returns {module:LinearAlgebra.SparseMatrix}
	 */
	static buildHodgeStar1Form(geometry, edgeIndex) {

		let E = geometry.mesh.edges.length;
		let T = new Triplet(E, E);
		let e, l_dual_over_l_primal;

		for (let i=0; i<E; i++) {
			e = geometry.mesh.edges[i];
			l_dual_over_l_primal = 0.5 * (geometry.cotan(e.halfedge) + geometry.cotan(e.halfedge.twin));
			T.addEntry(l_dual_over_l_primal, i, i);
		}

		let hodge1Form = SparseMatrix.fromTriplet(T);

		return hodge1Form;
	}

	/**
	 * Builds a sparse diagonal matrix encoding the Hodge operator on 2-forms.
	 * By convention, the area of a vertex is 1.
	 * @static
	 * @param {module:Core.Geometry} geometry The geometry of a mesh.
	 * @param {Object} faceIndex A dictionary mapping each face of a mesh to a unique index.
	 * @returns {module:LinearAlgebra.SparseMatrix}
	 */
	static buildHodgeStar2Form(geometry, faceIndex) {
		let F = geometry.mesh.faces.length;
		let T = new Triplet(F, F);
		let f, oneOverA;

		for (let i=0; i<F; i++) {
			f = geometry.mesh.faces[i];
			oneOverA = 1 / geometry.area(f);
			T.addEntry(oneOverA, i, i);
		}

		let hodge2Form = SparseMatrix.fromTriplet(T);

		return hodge2Form;
	}

	/**
	 * Builds a sparse matrix encoding the exterior derivative on 0-forms.
	 * @static
	 * @param {module:Core.Geometry} geometry The geometry of a mesh.
	 * @param {Object} edgeIndex A dictionary mapping each edge of a mesh to a unique index.
	 * @param {Object} vertexIndex A dictionary mapping each vertex of a mesh to a unique index.
	 * @returns {module:LinearAlgebra.SparseMatrix}
	 */
	static buildExteriorDerivative0Form(geometry, edgeIndex, vertexIndex) {
		let V = geometry.mesh.vertices.length;
		let E = geometry.mesh.edges.length;
		let T = new Triplet(E, V);
		let h, v_source, v_target;

		for (let i=0; i<E; i++) {
			h = geometry.mesh.edges[i].halfedge; 
			v_source = h.vertex.index;
			v_target = h.twin.vertex.index;
			T.addEntry(1, i, v_target);
			T.addEntry(-1, i, v_source);
		}

		return SparseMatrix.fromTriplet(T); // placeholder
	}

	/**
	 * Builds a sparse matrix encoding the exterior derivative on 1-forms.
	 * @static
	 * @param {module:Core.Geometry} geometry The geometry of a mesh.
	 * @param {Object} faceIndex A dictionary mapping each face of a mesh to a unique index.
	 * @param {Object} edgeIndex A dictionary mapping each edge of a mesh to a unique index.
	 * @returns {module:LinearAlgebra.SparseMatrix}
	 */
	static buildExteriorDerivative1Form(geometry, faceIndex, edgeIndex) {

		let E = geometry.mesh.edges.length;
		let F = geometry.mesh.faces.length;
		let T = new Triplet(F, E);
		let f, e, orientation;

		for (let i=0; i<F; i++) {

			f = geometry.mesh.faces[i];

			for (let h of f.adjacentHalfedges()) {

				e = geometry.mesh.edges[edgeIndex[h.edge.index]];
				if (h == e.halfedge) {  // check if current edge is the "primary" or twin edge of the halfedge
					orientation = 1;
				} else {
					orientation = -1;
				}

			T.addEntry(orientation, f.index, e.index);
			}
		}

		return SparseMatrix.fromTriplet(T);
	}
}
