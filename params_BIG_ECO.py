## For big single run on full distiller data set. load this file before running main.

organism = "Eco"

if ! isdefined(:k_clust) k_clust = 400; end

if ! isdefined(:ratios_file) ratios_file = "./junkey/$organism/DISTILLER_data.tsv"; end

if ! isdefined(:n_iters) n_iters = 200; end

## these are defaults (for microbes):
if ! isdefined(:distance_search) distance_search = [-150,+10]; end 
if ! isdefined(:distance_scan)   distance_scan =   [-200,+20]; end
if ! isdefined(:motif_width_range) motif_width_range = [6,30]; end

println( "$distance_search  $distance_scan  $motif_width_range" )

if ! isdefined(:max_network_weight) max_network_weight = 0.00001; end 
if ! isdefined(:max_motif_weight)   max_motif_weight =   0.8; end 

if ! isdefined(:avg_genes_per_cluster) avg_genes_per_cluster = 32; end 
if ! isdefined(:avg_clusters_per_gene) avg_clusters_per_gene = 12; end 

## allow more updates if there are more clusters??? Tuned to k_clust/2 for Hpy (where k_clust is 75) --
##     may need additional tuning; e.g. for eco (k_clust=450), k_clust/2 is too high
if ! isdefined(:max_improvements_per_iter) max_improvements_per_iter = convert(Int64, round(min(k_clust/2, 500))); end
