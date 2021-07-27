//----------------------------------------------------------------------
// File:			bd_search.cpp
// Programmer:		David Mount
// Description:		Standard bd-tree search
// Last modified:	01/04/05 (Version 1.0)
//----------------------------------------------------------------------
// Copyright (c) 1997-2005 University of Maryland and Sunil Arya and
// David Mount.  All Rights Reserved.
// 
// This software and related documentation is part of the Approximate
// Nearest Neighbor Library (ANN).  This software is provided under
// the provisions of the Lesser GNU Public License (LGPL).  See the
// file ../ReadMe.txt for further information.
// 
// The University of Maryland (U.M.) and the authors make no
// representations about the suitability or fitness of this software for
// any purpose.  It is provided "as is" without express or implied
// warranty.
//----------------------------------------------------------------------
// History:
//	Revision 0.1  03/04/98
//		Initial release
//----------------------------------------------------------------------

#include "bd_tree.h"					// bd-tree declarations
#include "kd_search.h"					// kd-tree search declarations

//----------------------------------------------------------------------
//	Approximate searching for bd-trees.
//		See the file kd_search.cpp for general information on the
//		approximate nearest neighbor search algorithm.  Here we
//		include the extensions for shrinking nodes.
//----------------------------------------------------------------------

//----------------------------------------------------------------------
//	bd_shrink::ann_search - search a shrinking node
//----------------------------------------------------------------------

void ANNbd_shrink::ann_search(ANNdist box_dist, int &ANNkdDim, ANNpoint &ANNkdQ, double &ANNkdMaxErr, ANNpointArray &ANNkdPts, ANNmin_k *ANNkdPointMK, int &ANNptsVisited)
{
												// check dist calc term cond.
	if (ANNmaxPtsVisited != 0 && ANNptsVisited > ANNmaxPtsVisited) return;

	ANNdist inner_dist = 0;						// distance to inner box
	for (int i = 0; i < n_bnds; i++) {			// is query point in the box?
		if (bnds[i].out(ANNkdQ)) {				// outside this bounding side?
												// add to inner distance
			inner_dist = (ANNdist) ANN_SUM(inner_dist, bnds[i].dist(ANNkdQ));
		}
	}
	if (inner_dist <= box_dist) {				// if inner box is closer
		child[ANN_IN]->ann_search(inner_dist, ANNkdDim, ANNkdQ, ANNkdMaxErr, ANNkdPts, ANNkdPointMK, ANNptsVisited);	// search inner child first
		child[ANN_OUT]->ann_search(box_dist, ANNkdDim, ANNkdQ, ANNkdMaxErr, ANNkdPts, ANNkdPointMK, ANNptsVisited);	// ...then outer child
	}
	else {										// if outer box is closer
		child[ANN_OUT]->ann_search(box_dist, ANNkdDim, ANNkdQ, ANNkdMaxErr, ANNkdPts, ANNkdPointMK, ANNptsVisited);	// search outer child first
		child[ANN_IN]->ann_search(inner_dist, ANNkdDim, ANNkdQ, ANNkdMaxErr, ANNkdPts, ANNkdPointMK, ANNptsVisited);	// ...then outer child
	}
	ANN_FLOP(3*n_bnds)							// increment floating ops
	ANN_SHR(1)									// one more shrinking node
}

//----------------------------------------------------------------------
//	annkSearch - search for the k nearest neighbors
//----------------------------------------------------------------------

void ANNbd_tree::annkSearch(
	ANNpoint			q,				// the query point
	int					k,				// number of near neighbors to return
	ANNidxArray			nn_idx,			// nearest neighbor indices (returned)
	ANNdistArray		dd,				// the approximate nearest neighbor
	double				eps)			// the error bound
{

    int				ANNkdDim = dim;		// dimension of space (static copy)
    ANNpoint	    ANNkdQ = q;			// query point (static copy)
    double			ANNkdMaxErr;	    // max tolerable squared error
    ANNpointArray	ANNkdPts = pts;		// the points (static copy)
    ANNmin_k	    *ANNkdPointMK;	    // set of k closest points
    int				ANNptsVisited = 0;	// number of points visited


	if (k > n_pts) {					// too many near neighbors?
		annError("Requesting more near neighbors than data points", ANNabort);
	}

	ANNkdMaxErr = ANN_POW(1.0 + eps);
	ANN_FLOP(2)							// increment floating op count

	ANNkdPointMK = new ANNmin_k(k);		// create set for closest k points
										// search starting at the root
	root->ann_search(annBoxDistance(q, bnd_box_lo, bnd_box_hi, dim), ANNkdDim, ANNkdQ, ANNkdMaxErr, ANNkdPts, ANNkdPointMK, ANNptsVisited);

	for (int i = 0; i < k; i++) {		// extract the k-th closest points
		dd[i] = ANNkdPointMK->ith_smallest_key(i);
		nn_idx[i] = ANNkdPointMK->ith_smallest_info(i);
	}
	delete ANNkdPointMK;				// deallocate closest point set
}
