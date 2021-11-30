library("tidyverse")
library("abcrf")
library("ggthemes")
library("scales")

submodels = c("AM_1M_1N", "AM_1M_2N", "AM_2M_1N", "AM_2M_2N", "IM_1M_1N", "IM_1M_2N", "IM_2M_1N", "IM_2M_2N", "SC_1M_1N", "SC_1M_2N", "SC_2M_1N", "SC_2M_2N", "SI_1N", "SI_2N")
models = c("isolation", "isolation", "isolation", "isolation", "migration", "migration", "migration", "migration", "migration", "migration", "migration", "migration", "isolation", "isolation")

# 1: make the reference
sim_ref = NULL
index_submodels = c() # in order to know which simulation/iteration comes from which submodel
index_models = c()
for(i in 1:length(submodels)){
	submodel_tmp = submodels[i]
	for(iteration in dir(path="reference", pattern=submodel_tmp)){
		# get the stats
		tmp_stats = read.table(paste("reference", iteration, "ABCstat.txt", sep="/"), h=T)
		
		# binds the tables
		sim_ref = rbind(sim_ref, tmp_stats)
		
		# fill the index vectors
		index_submodels = c(index_submodels, rep(x=submodel_tmp, times=nrow(tmp_stats)))
		index_models = c(index_models, rep(x=models[i], times=nrow(tmp_stats)))
	}
}

# 2: grow the forest
index_models=as.factor(index_models)
data1 = data.frame(index_models, sim_ref[,-(1:3)]) # remove the first column because it's useless, and col2+3 because they are invariable
model_rf = abcrf(index_models~., data = data1, ntree=1000)

# 3: get the PODs (for Pseudo Observed datasets)
sim_pods = NULL
index_submodels_pods = c() # in order to know which simulation/iteration comes from which submodel
index_models_pods = c()
for(i in 1:length(submodels)){
	submodel_tmp = submodels[i]
	for(iteration in dir(path="PODs", pattern=submodel_tmp)){
		# get the stats
		tmp_stats = read.table(paste("PODs", iteration, "ABCstat.txt", sep="/"), h=T)
		
		# binds the tables
		sim_pods = rbind(sim_pods, tmp_stats)
		
		# fill the index vectors
		index_submodels_pods = c(index_submodels_pods, rep(x=submodel_tmp, times=nrow(tmp_stats)))
		index_models_pods = c(index_models_pods, rep(x=models[i], times=nrow(tmp_stats)))
	}
}

# 4: apply the forest on the PODs
prediction_PODs = predict(model_rf, sim_pods[, -(1:3)], data1, ntree=1000)

# 5: be an artist
couleurs = c("#3B9AB2", "#9EBE91", "#E4B80E", "#F21A00")
res = tibble(submodel=index_submodels_pods, real_model=index_models_pods, estimated_model=prediction_PODs$allocation, post_proba_bestModel=prediction_PODs$post.prob, netDivergence=sim_pods$netdivAB_avg)
res = res %>% dplyr::mutate(proba_migration=ifelse(estimated_model=="migration", post_proba_bestModel, 1-post_proba_bestModel))
res = res %>% dplyr::mutate(submodel_demo=factor(substr(x=submodel, start=1, stop=2), levels=c('SI', 'AM', 'SC', 'IM')))

P_grey_zone = res %>% ggplot(aes(x=netDivergence, y=proba_migration, col=submodel_demo)) +     
	geom_point(size=3, alpha=0.15) + facet_grid(cols=vars(submodel_demo)) +  
	xlim(1e-5, 0.3) + ylim(0, 1) + theme_calc(base_size = 17) + 
	labs(x="net divergence", y=expression(P["migration"])) + 
	scale_color_manual(values=couleurs, "Model") +
	scale_x_log10( breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) +
	guides(colour=guide_legend(override.aes=list(alpha=1)))

ggsave(filename="plot_proba_models.pdf", plot=P_grey_zone, device="pdf", width=17, height=4.5, units="in", bg="white")

# 6 quantifying the observed plot
## explore the continuum of divergence
minDiv_list = c(0, 1e-4, 5e-4, 1e-3, 5e-3, 1e-2, 5e-2, 1e-1)
maxDiv_list = c(1e-4, 5e-4, 1e-3, 5e-3, 1e-2, 5e-2, 1e-1, 1)
res_cnt = NULL
levels_bin = NULL
for(i in 1:length(minDiv_list)){
	minDiv = minDiv_list[i]
	maxDiv = maxDiv_list[i]
	bin_name=paste(minDiv, "<Da<", maxDiv, sep="")
	levels_bin = c(levels_bin, bin_name)
	tmp = res %>% dplyr::filter(netDivergence>minDiv, netDivergence<=maxDiv) %>% group_by(submodel_demo, estimated_model, .drop=F) %>% summarise(count=n())
	tmp = tmp %>% dplyr::mutate(bin=bin_name)
	res_cnt = rbind(res_cnt, tmp)
}

res_cnt = res_cnt %>% dplyr::mutate(bin=factor(bin, levels=levels_bin))

P_nSim_divergence = res_cnt %>% group_by(submodel_demo, estimated_model) %>% ggplot(aes(x=bin, y=count, fill=submodel_demo)) + geom_bar(stat="identity", position=position_dodge()) +
	theme_calc(base_size = 17) + 
	labs(x="net divergence", y="count") + 
	scale_fill_manual(values=couleurs, "Model")

ggsave(filename="plot_nSim.pdf", plot=P_nSim_divergence, device="pdf", width=17, height=6, units="in", bg="white")

