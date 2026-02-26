"use strict";

/**
 * @module Projects
 */
class SimplicialComplexOperators {

        /** This class implements various operators (e.g. boundary, star, link) on a mesh.
         * @constructor module:Projects.SimplicialComplexOperators
         * @param {module:Core.Mesh} mesh The input mesh this class acts on.
         * @property {module:Core.Mesh} mesh The input mesh this class acts on.
         * @property {module:LinearAlgebra.SparseMatrix} A0 The vertex-edge adjacency matrix of <code>mesh</code>.
         * @property {module:LinearAlgebra.SparseMatrix} A1 The edge-face adjacency matrix of <code>mesh</code>.
         */
        constructor(mesh) {
                this.mesh = mesh;
                this.assignElementIndices(this.mesh);

                this.A0 = this.buildVertexEdgeAdjacencyMatrix(this.mesh);
                this.A1 = this.buildEdgeFaceAdjacencyMatrix(this.mesh);
        }

        /** Assigns indices to the input mesh's vertices, edges, and faces
         * @method module:Projects.SimplicialComplexOperators#assignElementIndices
         * @param {module:Core.Mesh} mesh The input mesh which we index.
         */
        assignElementIndices(mesh) {
                for (let i=0; i<mesh.vertices.length; i++) {
                        let v = mesh.vertices[i];
                        v.index = i;
                }
                for (let i=0; i<mesh.edges.length; i++) {
                        mesh.edges[i].index = i;
                }
                for (let i=0; i<mesh.faces.length; i++) {
                        mesh.faces[i].index = i;
                }
        }

        /** Returns the vertex-edge adjacency matrix of the given mesh.
         * @method module:Projects.SimplicialComplexOperators#buildVertexEdgeAdjacencyMatrix
         * @param {module:Core.Mesh} mesh The mesh whose adjacency matrix we compute.
         * @returns {module:LinearAlgebra.SparseMatrix} The vertex-edge adjacency matrix of the given mesh.
         */
        buildVertexEdgeAdjacencyMatrix(mesh) {

                let n_v = mesh.vertices.length;
                let n_e = mesh.edges.length;
                let t_edges = new Triplet(n_e, n_v);
                
                for (let i=0; i<n_e; i++) {
                        let e = mesh.edges[i];
                        let e_ind = e.index;
                        let v_source = e.halfedge.vertex.index;
                        let v_dest = e.halfedge.twin.vertex.index;
                        t_edges.addEntry(1, e_ind, v_source);
                        t_edges.addEntry(1, e_ind, v_dest);
                }

                let A0 = SparseMatrix.fromTriplet(t_edges);

                return A0;
        }
                
        /** Returns the edge-face adjacency matrix.
         * @param {module:Core.Mesh} mesh The mesh whose adjacency matrix we compute.
         * @returns {module:LinearAlgebra.SparseMatrix} The edge-face adjacency matrix of the given mesh.
         */
        buildEdgeFaceAdjacencyMatrix(mesh) {

                let n_e = mesh.edges.length;
                let n_f = mesh.faces.length;
                let t_faces = new Triplet(n_f, n_e);
                
                for (let i=0; i<n_f; i++) {
                        let f = mesh.faces[i];
                        let f_ind = f.index;
                        let e_1 = f.halfedge.edge.index;
                        let e_2 = f.halfedge.next.edge.index;
                        let e_3 = f.halfedge.next.next.edge.index;
                        t_faces.addEntry(1, f_ind, e_1);
                        t_faces.addEntry(1, f_ind, e_2);
                        t_faces.addEntry(1, f_ind, e_3);
                }

                let A1 = SparseMatrix.fromTriplet(t_faces);

                return A1;
        }

        /** Returns a column vector representing the vertices of the
         * given subset.
         * @method module:Projects.SimplicialComplexOperators#buildVertexVector
         * @param {module:Core.MeshSubset} subset A subset of our mesh.
         * @returns {module:LinearAlgebra.DenseMatrix} A column vector with |V| entries. The ith entry is 1 if
         *  vertex i is in the given subset and 0 otherwise
         */
        buildVertexVector(subset) {

                let V = this.mesh.vertices.length;
                let vec = DenseMatrix.zeros(V, 1);
                for (let val of subset.vertices) {
                        vec.set(1, val, 0);
                }

                return vec;
        }

        /** Returns a column vector representing the edges of the
         * given subset.
         * @method module:Projects.SimplicialComplexOperators#buildEdgeVector
         * @param {module:Core.MeshSubset} subset A subset of our mesh.
         * @returns {module:LinearAlgebra.DenseMatrix} A column vector with |E| entries. The ith entry is 1 if
         *  edge i is in the given subset and 0 otherwise
         */
        buildEdgeVector(subset) {
                
                let E = this.mesh.edges.length;
                let vec = DenseMatrix.zeros(E, 1);
                for (let val of subset.edges) {
                        vec.set(1, val, 0);
                }

                return vec;
        }

        /** Returns a column vector representing the faces of the
         * given subset.
         * @method module:Projects.SimplicialComplexOperators#buildFaceVector
         * @param {module:Core.MeshSubset} subset A subset of our mesh.
         * @returns {module:LinearAlgebra.DenseMatrix} A column vector with |F| entries. The ith entry is 1 if
         *  face i is in the given subset and 0 otherwise
         */
        buildFaceVector(subset) {

                let F = this.mesh.faces.length;
                let vec = DenseMatrix.zeros(F, 1);
                for (let val of subset.faces) {
                        vec.set(1, val, 0);
                }

                return vec;
        }

        /** Returns the star of a subset.
         * @method module:Projects.SimplicialComplexOperators#star
         * @param {module:Core.MeshSubset} subset A subset of our mesh.
         * @returns {module:Core.MeshSubset} The star of the given subset.
         */
        star(subset) {
                /*
                Star of a mesh contains its own simplices (i.e. all vertices, edges, faces) and 
                all of the simplices that contain any of these simplices.
                
                We can use the adjaceny metrices in order to get all the needed simplices.
                        vertices_vector will give all the input vertices.
                        edges_vector + A0 @ vertices_vector will give all the edges that contain any of the vertices
                        faces_vector + A1 @ edges_vector + A1 @ A0 @ vertices_vector will give all the faces that contain any of the vertices or edges

                reference - Gemini answer:
                https://gemini.google.com/share/f0fcc821c5ef
                */
                
                let vertices_vector = this.buildVertexVector(subset);
                let edges_vector = this.buildEdgeVector(subset);
                let faces_vector = this.buildFaceVector(subset);

                let v_star = vertices_vector;
                let e_star = edges_vector.plus(this.A0.timesDense(vertices_vector));
                let f_star = faces_vector.plus(this.A1.timesDense(e_star));

                let v_star_set = new Set();
                for (let i=0; i<vertices_vector.nRows(); i++) {
                        if (v_star.get(i, 0) > 0) {
                                v_star_set.add(i);
                        }
                }
                
                let e_star_set = new Set();
                for (let i=0; i<edges_vector.nRows(); i++) {
                        if (e_star.get(i, 0) > 0) {
                                e_star_set.add(i);
                        }
                }
                
                let f_star_set = new Set();
                for (let i=0; i<faces_vector.nRows(); i++) {
                        if (f_star.get(i, 0) > 0) {
                                f_star_set.add(i);
                        }
                }

                let star = new MeshSubset(v_star_set, e_star_set, f_star_set);
                
                
                return star;
        }

        /** Returns the closure of a subset.
         * @method module:Projects.SimplicialComplexOperators#closure
         * @param {module:Core.MeshSubset} subset A subset of our mesh.
         * @returns {module:Core.MeshSubset} The closure of the given subset.
         */
        closure(subset) {
                /*
                function logic:
                        - loop over all faces and add all edges 
                        - loop over all edges (include those who were just added) and add all vertices
                */

                let edges_vector = this.buildEdgeVector(subset);
                let faces_vector = this.buildFaceVector(subset);
             
                let e_closure = edges_vector.plus(this.A1.transpose().timesDense(faces_vector));
                
                let e_closure_set = new Set();
                for (let i=0; i<e_closure.nRows(); i++) {
                        if (e_closure.get(i, 0) > 0) {
                                e_closure_set.add(i);
                        }
                }

                let vertices_vector = this.buildVertexVector(subset);
                let v_closure = vertices_vector.plus(this.A0.transpose().timesDense(e_closure));

                let v_closure_set = new Set();
                for (let i=0; i<v_closure.nRows(); i++) {
                        if (v_closure.get(i, 0) > 0) {
                                v_closure_set.add(i);
                        }
                }
                
                let closure = new MeshSubset(v_closure_set, e_closure_set, subset.faces);
                
                return closure;
        }

        /** Returns the link of a subset.
         * @method module:Projects.SimplicialComplexOperators#link
         * @param {module:Core.MeshSubset} subset A subset of our mesh.
         * @returns {module:Core.MeshSubset} The link of the given subset.
         */
        link(subset) {

                let link = new MeshSubset
                link = this.closure(this.star(subset))
                link.deleteSubset(this.star(this.closure(subset)))  // inplace

                return link;
        }

        /** Returns true if the given subset is a subcomplex and false otherwise.
         * @method module:Projects.SimplicialComplexOperators#isComplex
         * @param {module:Core.MeshSubset} subset A subset of our mesh.
         * @returns {boolean} True if the given subset is a subcomplex and false otherwise.
         */
        isComplex(subset) {

                let subset_copy = new MeshSubset(new Set(subset.vertices), new Set(subset.edges), new Set(subset.faces));  // make a copy since deleteSubset() deletes faces from input subset
                let closure = new MeshSubset
                closure = this.closure(subset_copy)
                closure.deleteSubset(subset_copy)  // inplace

                let is_complex = (closure.vertices.size == 0) && (closure.edges.size == 0) && (closure.faces.size == 0);
                return is_complex;
        }

        /** Returns the degree if the given subset is a pure subcomplex and -1 otherwise.
         * @method module:Projects.SimplicialComplexOperators#isPureComplex
         * @param {module:C. re.MeshSubset} subset A subset of our mesh.
         * @returns {number} The degree of the given subset if it is a pure subcomplex and -1 otherwise.
         */
        isPureComplex(subset) {
                /* chapter 2.1 from the notes
                A complex K is a pure k-simplicial complex if every simplex σ ∈ K is contained in some simplex of degree k (possibly itself). 
                For instance, a bunch of triangles with edges and vertices hanging off the side or floating around by themselves is not pure.
                */

                if (!this.isComplex(subset)) {
                        return -1;
                }


                if (subset.faces.size > 0) {
                        
                        let faces_vector = this.buildFaceVector(subset);
                        let vertices_of_faces = this.A0.transpose().timesDense(this.A1.transpose().timesDense(faces_vector))
                        let vertices_vector = this.buildVertexVector(subset)
                     // compare vectors
                        for (let i=0; i<vertices_of_faces.nRows(); i++) {
                                let flag = ((vertices_of_faces.get(i, 0) > 0) && (vertices_vector.get(i, 0) > 0)) || ((vertices_of_faces.get(i, 0) == 0) && (vertices_vector.get(i, 0) == 0));
                                if (!flag) {
                                        return -1;
                                }
                        }
                        
                        return 2;
           
                        
                } else if (subset.edges.size > 0) {
                        
                        let edges_vector = this.buildEdgeVector(subset);
                        let vertices_of_edges = this.A0.transpose().timesDense(edges_vector)
                        let vertices_vector = this.buildVertexVector(subset)
                        // compare vectors
                        for (let i=0; i<vertices_of_edges.nRows(); i++) {
                                let flag = ((vertices_of_edges.get(i, 0) > 0) && (vertices_vector.get(i, 0) > 0)) || ((vertices_of_edges.get(i, 0) == 0) && (vertices_vector.get(i, 0) == 0));
                                if (!flag) {
                                        return -1;
                                }
                        }
                        
                        return 1;
                        
                } else if (subset.vertices.size > 0) { 
                        // only edges - they are always pure complex
                        return 0;
                }
                
                return -1;
        }

        /** Returns the boundary of a subset.
         * @method module:Projects.SimplicialComplexOperators#boundary
         * @param {module:Core.MeshSubset} subset A subset of our mesh. We assume <code>subset</code> is a pure subcomplex.
         * @returns {module:Core.MeshSubset} The boundary of the given pure subcomplex.
         */
        boundary(subset) {
                // TODO

                return subset; // placeholder
        }
}
