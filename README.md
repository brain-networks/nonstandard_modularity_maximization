# nonstandard_modularity_maximization

This repository contains a series of scripts for generating non-standard null models to be used with modularity maximization.

What is here?
1. single_layer_null_models.m - generates null models and modularity matrices for single-layer networks (see Fig. 2b and Fig. 2c).

  - example 1. Newman-Girvan degree-/strength-preserving null model (with self-connections)
	
  - example 2. Newman-Girvan degree-/strength-preserving null model (without self-connections)
	
  - example 3. uniform null model where all connections (except for self-loops) have an expected weight equal to that of the resolution parameter, ``gamma''
	
  - example 4. geometric (spatial) model where connections are drawn stochastically from a decaying exponential. weights are added to the connections to be inversely proportional to distance.
	
  - example 5. minimally wired network containing the shortest possible set of connections. weights are added to the connections to be inversely proportional to distance.
	
  - example 6. spatial model that also preserves binary degree sequence (but not strength).
	
  - example 7. preserves the binary topology of the network but assigns edges a uniform weight equal to the mean weight across all edges.
	
  - example 8. signed and weighted version of the Newman-Girvan null model. this version weights the contribution of positive/negative connections equally.
	
  - example 9. signed and weighted version of the Newman-Girvan null model. this version weights positive contributions more strongly.
		
2. multi_layer_examples.m - generates flattened multi-layer modularity matrices (see Fig. 4).
3. fcn/ - a set of helper functions for implementing the different models.
4. data/ - example FC and SC matrices.

If you use any of these scripts, please cite:
Esfahlani, F. Z., Jo, Y., Puxeddu, M. G., Merritt, H., Tanner, J. C., Greenwell, S., ... & Betzel, R. F. (2021). Modularity maximization as a flexible and generic framework for brain network exploratory analysis. arXiv preprint arXiv:2106.15428.