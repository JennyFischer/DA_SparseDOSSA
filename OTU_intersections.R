#relab intersection doesnt work properly and gets a list of zero intersections 
#altough there are some same OTUs
#it just works when doing an intersection with another list of relab
#wilcoxon
#clr
wlx_lognormal_norm_clr<- read.csv("/Users/jf/Desktop/BA/Analysis_SparseDOSSA/outputs/wilcoxon_output/wlx_pvalues_lognormal_norm_clr.csv")
wlx_lognormal_norm_clr_outs = as.list(subset(wlx_lognormal_norm_clr[wlx_lognormal_norm_clr[,"p"] < 0.05, ], wlx_lognormal_norm_clr[wlx_lognormal_norm_clr[,"p"] < 0.05, ][1] != "NA")[1])
wlx_lognormal_norm_clr_outs
#css
wlx_lognormal_norm_css<- read.csv("/Users/jf/Desktop/BA/Analysis_SparseDOSSA/outputs/wilcoxon_output/wlx_pvalues_lognormal_norm_css.csv")
wlx_lognormal_norm_css_outs = as.list(subset(wlx_lognormal_norm_css[wlx_lognormal_norm_css[,"p"] < 0.05, ], wlx_lognormal_norm_css[wlx_lognormal_norm_css[,"p"] < 0.05, ][1] != "NA")[1])
wlx_lognormal_norm_css_outs
#tss
wlx_lognormal_norm_tss<- read.csv("/Users/jf/Desktop/BA/Analysis_SparseDOSSA/outputs/wilcoxon_output/wlx_pvalues_lognormal_norm_tss.csv")
wlx_lognormal_norm_tss_outs = as.list(subset(wlx_lognormal_norm_tss[wlx_lognormal_norm_tss[,"p"] < 0.05, ], wlx_lognormal_norm_tss[wlx_lognormal_norm_tss[,"p"] < 0.05, ][1] != "NA")[1])
wlx_lognormal_norm_tss_outs
#relab
wlx_lognormal_norm_relab<- read.csv("/Users/jf/Desktop/BA/Analysis_SparseDOSSA/outputs/wilcoxon_output/wlx_pvalues_lognormal_norm_relab.csv")
wlx_lognormal_norm_relab_outs = as.list(subset(wlx_lognormal_norm_relab[wlx_lognormal_norm_relab[,"p"] < 0.05, ], wlx_lognormal_norm_relab[wlx_lognormal_norm_relab[,"p"] < 0.05, ][1] != "NA")[1])
wlx_lognormal_norm_relab_outs

#spike
#clr
wlx_logspike_norm_clr<- read.csv("/Users/jf/Desktop/BA/Analysis_SparseDOSSA/outputs/wilcoxon_output/wlx_pvalues_logspike_norm_clr.csv")
wlx_logspike_norm_clr_outs = as.list(subset(wlx_logspike_norm_clr[wlx_logspike_norm_clr[,"p"] < 0.05, ], wlx_logspike_norm_clr[wlx_logspike_norm_clr[,"p"] < 0.05, ][1] != "NA")[1])
wlx_logspike_norm_clr_outs
#css
wlx_logspike_norm_css<- read.csv("/Users/jf/Desktop/BA/Analysis_SparseDOSSA/outputs/wilcoxon_output/wlx_pvalues_logspike_norm_css.csv")
wlx_logspike_norm_css_outs = as.list(subset(wlx_logspike_norm_css[wlx_logspike_norm_css[,"p"] < 0.05, ], wlx_logspike_norm_css[wlx_logspike_norm_css[,"p"] < 0.05, ][1] != "NA")[1])
wlx_logspike_norm_css_outs
#tss
wlx_logspike_norm_tss<- read.csv("/Users/jf/Desktop/BA/Analysis_SparseDOSSA/outputs/wilcoxon_output/wlx_pvalues_logspike_norm_tss.csv")
wlx_logspike_norm_tss_outs = as.list(subset(wlx_logspike_norm_tss[wlx_logspike_norm_tss[,"p"] < 0.05, ], wlx_logspike_norm_tss[wlx_logspike_norm_tss[,"p"] < 0.05, ][1] != "NA")[1])
wlx_logspike_norm_tss_outs
#relab
wlx_logspike_norm_relab<- read.csv("/Users/jf/Desktop/BA/Analysis_SparseDOSSA/outputs/wilcoxon_output/wlx_pvalues_logspike_norm_relab.csv")
wlx_logspike_norm_relab_outs = as.list(subset(wlx_logspike_norm_relab[wlx_logspike_norm_relab[,"p"] < 0.05, ], wlx_logspike_norm_relab[wlx_logspike_norm_relab[,"p"] < 0.05, ][1] != "NA")[1])
wlx_logspike_norm_relab_outs




intersect_wlx_lognormal_norm_css_clr = intersect(wlx_lognormal_norm_clr_outs,wlx_lognormal_norm_css_outs)
intersect_wlx_lognormal_norm_tss_clr = intersect(wlx_lognormal_norm_clr_outs,wlx_lognormal_norm_tss_outs)
intersect_wlx_lognormal_norm_css_tss = intersect(wlx_lognormal_norm_css_outs,wlx_lognormal_norm_tss_outs)
intersect_wlx_lognormal_norm_rel_clr = intersect(wlx_lognormal_norm_clr_outs,wlx_lognormal_norm_relab_outs)
intersect_wlx_lognormal_norm_rel_tss = intersect(wlx_lognormal_norm_relab_outs,wlx_lognormal_norm_tss_outs)
intersect_wlx_lognormal_norm_css_rel = intersect(wlx_lognormal_norm_css_outs,wlx_lognormal_norm_relab_outs)
intersect_wlx_lognormal_norm_css_clr
intersect_wlx_lognormal_norm_tss_clr

intersect_wlx_logspike_norm_css_clr = intersect(wlx_logspike_norm_clr_outs,wlx_logspike_norm_css_outs)
intersect_wlx_logspike_norm_tss_clr = intersect(wlx_logspike_norm_clr_outs,wlx_logspike_norm_tss_outs)
intersect_wlx_logspike_norm_css_tss = intersect(wlx_logspike_norm_css_outs,wlx_logspike_norm_tss_outs)
intersect_wlx_logspike_norm_rel_clr = intersect(wlx_logspike_norm_clr_outs,wlx_logspike_norm_relab_outs)
intersect_wlx_logspike_norm_rel_tss = intersect(wlx_logspike_norm_relab_outs,wlx_logspike_norm_tss_outs)
intersect_wlx_logspike_norm_css_rel = intersect(wlx_logspike_norm_css_outs,wlx_logspike_norm_relab_outs)
intersect_wlx_logspike_norm_css_clr
intersect_wlx_logspike_norm_tss_clr
intersect_wlx_logspike_norm_css_tss
intersect_wlx_logspike_norm_rel_clr #doesnt work correctly
intersect_wlx_logspike_norm_rel_tss #doesnt work correctly
intersect_wlx_logspike_norm_css_rel #doesnt work correctly

#Ttest
#clr
ttest_lognormal_norm_clr<- read.csv("/Users/jf/Desktop/BA/Analysis_SparseDOSSA/outputs/ttest_output/tt_pvalues_lognormal_norm_clr.csv")
ttest_lognormal_norm_clr_outs = as.list(subset(wlx_lognormal_norm_clr[wlx_lognormal_norm_clr[,"p"] < 0.05, ], wlx_lognormal_norm_clr[wlx_lognormal_norm_clr[,"p"] < 0.05, ][1] != "NA")[1])
ttest_lognormal_norm_clr_outs
#css
ttest_lognormal_norm_css<- read.csv("/Users/jf/Desktop/BA/Analysis_SparseDOSSA/outputs/ttest_output/tt_pvalues_lognormal_norm_css.csv")
ttest_lognormal_norm_css_outs = as.list(subset(wlx_lognormal_norm_css[wlx_lognormal_norm_css[,"p"] < 0.05, ], wlx_lognormal_norm_css[wlx_lognormal_norm_css[,"p"] < 0.05, ][1] != "NA")[1])
ttest_lognormal_norm_css_outs
#tss
ttest_lognormal_norm_tss<- read.csv("/Users/jf/Desktop/BA/Analysis_SparseDOSSA/outputs/ttest_output/tt_pvalues_lognormal_norm_tss.csv")
ttest_lognormal_norm_tss_outs = as.list(subset(wlx_lognormal_norm_tss[wlx_lognormal_norm_tss[,"p"] < 0.05, ], wlx_lognormal_norm_tss[wlx_lognormal_norm_tss[,"p"] < 0.05, ][1] != "NA")[1])
ttest_lognormal_norm_tss_outs
#relab
ttest_lognormal_norm_relab<- read.csv("/Users/jf/Desktop/BA/Analysis_SparseDOSSA/outputs/ttest_output/tt_pvalues_lognormal_norm_relab.csv")
ttest_lognormal_norm_relab_outs = as.list(subset(wlx_lognormal_norm_relab[wlx_lognormal_norm_relab[,"p"] < 0.05, ], wlx_lognormal_norm_relab[wlx_lognormal_norm_relab[,"p"] < 0.05, ][1] != "NA")[1])
ttest_lognormal_norm_relab_outs

#ttest spike
#clr
ttest_logspike_norm_clr<- read.csv("/Users/jf/Desktop/BA/Analysis_SparseDOSSA/outputs/ttest_output/tt_pvalues_logspike_norm_clr.csv")
ttest_logspike_norm_clr_outs = as.list(subset(wlx_logspike_norm_clr[wlx_logspike_norm_clr[,"p"] < 0.05, ], wlx_logspike_norm_clr[wlx_logspike_norm_clr[,"p"] < 0.05, ][1] != "NA")[1])
ttest_logspike_norm_clr_outs
#css
ttest_logspike_norm_css<- read.csv("/Users/jf/Desktop/BA/Analysis_SparseDOSSA/outputs/ttest_output/tt_pvalues_logspike_norm_css.csv")
ttest_logspike_norm_css_outs = as.list(subset(wlx_logspike_norm_css[wlx_logspike_norm_css[,"p"] < 0.05, ], wlx_logspike_norm_css[wlx_logspike_norm_css[,"p"] < 0.05, ][1] != "NA")[1])
ttest_logspike_norm_css_outs
#tss
ttest_logspike_norm_tss<- read.csv("/Users/jf/Desktop/BA/Analysis_SparseDOSSA/outputs/ttest_output/tt_pvalues_logspike_norm_tss.csv")
ttest_logspike_norm_tss_outs = as.list(subset(wlx_logspike_norm_tss[wlx_logspike_norm_tss[,"p"] < 0.05, ], wlx_logspike_norm_tss[wlx_logspike_norm_tss[,"p"] < 0.05, ][1] != "NA")[1])
ttest_logspike_norm_tss_outs
#relab
ttest_logspike_norm_relab<- read.csv("/Users/jf/Desktop/BA/Analysis_SparseDOSSA/outputs/ttest_output/tt_pvalues_logspike_norm_relab.csv")
ttest_logspike_norm_relab_outs = as.list(subset(wlx_logspike_norm_relab[wlx_logspike_norm_relab[,"p"] < 0.05, ], wlx_logspike_norm_relab[wlx_logspike_norm_relab[,"p"] < 0.05, ][1] != "NA")[1])
ttest_logspike_norm_relab_outs 


intersect_ttest_lognormal_norm_css_clr = intersect(ttest_lognormal_norm_clr_outs,ttest_lognormal_norm_css_outs)
intersect_ttest_lognormal_norm_tss_clr = intersect(ttest_lognormal_norm_clr_outs,ttest_lognormal_norm_tss_outs)
intersect_ttest_lognormal_norm_css_tss = intersect(ttest_lognormal_norm_css_outs,ttest_lognormal_norm_tss_outs)
intersect_ttest_lognormal_norm_rel_clr = intersect(ttest_lognormal_norm_clr_outs,ttest_lognormal_norm_relab_outs)
intersect_ttest_lognormal_norm_rel_tss = intersect(ttest_lognormal_norm_relab_outs,ttest_lognormal_norm_tss_outs)
intersect_ttest_lognormal_norm_css_rel = intersect(ttest_lognormal_norm_css_outs,ttest_lognormal_norm_relab_outs)
intersect_ttest_lognormal_norm_css_clr
intersect_ttest_lognormal_norm_tss_clr
intersect_ttest_lognormal_norm_css_tss
intersect_ttest_lognormal_norm_rel_clr
intersect_ttest_lognormal_norm_rel_tss
intersect_ttest_lognormal_norm_css_rel

intersect_ttest_logspike_norm_css_clr = intersect(ttest_logspike_norm_clr_outs,ttest_logspike_norm_css_outs)
intersect_ttest_logspike_norm_tss_clr = intersect(ttest_logspike_norm_clr_outs,ttest_logspike_norm_tss_outs)
intersect_ttest_logspike_norm_css_tss = intersect(ttest_logspike_norm_css_outs,ttest_logspike_norm_tss_outs)
intersect_ttest_logspike_norm_rel_clr = intersect(ttest_logspike_norm_clr_outs,ttest_logspike_norm_relab_outs)
intersect_ttest_logspike_norm_rel_tss = intersect(ttest_logspike_norm_relab_outs,ttest_logspike_norm_tss_outs)
intersect_ttest_logspike_norm_css_rel = intersect(ttest_logspike_norm_css_outs,ttest_logspike_norm_relab_outs)

intersect_ttest_logspike_norm_css_clr
intersect_ttest_logspike_norm_tss_clr
intersect_ttest_logspike_norm_css_tss
intersect_ttest_logspike_norm_rel_clr
intersect_ttest_logspike_norm_rel_tss
intersect_ttest_logspike_norm_css_rel 


#aldex2
#lognorm welchs test
aldex2_lognormal_wel<- read.csv("/Users/jf/Desktop/BA/Analysis_SparseDOSSA/outputs/aldex2_output/aldex2_lognormal_welch.csv")
aldex2_lognormal_wel_outs = as.list(subset(aldex2_lognormal_wel[aldex2_lognormal_wel[,"we.ep"] < 0.05, ], aldex2_lognormal_wel[aldex2_lognormal_wel[,"we.ep"] < 0.05, ][1] != "NA")[2])
aldex2_lognormal_wel_outs
#lognorm wilcoxon
aldex2_lognormal_wil<- read.csv("/Users/jf/Desktop/BA/Analysis_SparseDOSSA/outputs/aldex2_output/aldex2_lognormal_wil.csv")
aldex2_lognormal_wil_outs = as.list(subset(aldex2_lognormal_wil[aldex2_lognormal_wil[,"wi.ep"] < 0.05, ], aldex2_lognormal_wil[aldex2_lognormal_wil[,"wi.ep"] < 0.05, ][1] != "NA")[2])
aldex2_lognormal_wil_outs
#spike welchs test
aldex2_spike_wel<- read.csv("/Users/jf/Desktop/BA/Analysis_SparseDOSSA/outputs/aldex2_output/aldex2_logspike_welch.csv")
aldex2_spike_wel_outs = as.list(subset(aldex2_spike_wel[aldex2_spike_wel[,"we.ep"] < 0.05, ], aldex2_spike_wel[aldex2_spike_wel[,"we.ep"] < 0.05, ][1] != "NA")[2])
aldex2_spike_wel_outs
#spike wilcoxon
aldex2_spike_wil<- read.csv("/Users/jf/Desktop/BA/Analysis_SparseDOSSA/outputs/aldex2_output/aldex2_logspike_wil.csv")
aldex2_spike_wil_outs = as.list(subset(aldex2_spike_wil[aldex2_spike_wil[,"wi.ep"] < 0.05, ], aldex2_spike_wil[aldex2_spike_wil[,"wi.ep"] < 0.05, ][1] != "NA")[2])
aldex2_spike_wil_outs


intersect_aldex2_wel = intersect(aldex2_lognormal_wel_outs,aldex2_spike_wel_outs)
intersect_aldex2_wil = intersect(aldex2_lognormal_wil_outs,aldex2_spike_wil_outs)
intersect_aldex2_wel = intersect(aldex2_lognormal_wel_outs,aldex2_lognormal_wil_outs)
intersect_aldex2_wel = intersect(aldex2_spike_wel_outs,aldex2_spike_wil_outs)

#Wilcoxon
#normal
#clr
wlx_lognormal_norm_clr<- read.csv("/Users/jf/Desktop/BA/Analysis_SparseDOSSA/outputs/wilcoxon_output/wlx_pvalues_lognormal_norm_clr.csv")
wlx_lognormal_norm_clr_outs = as.list(subset(wlx_lognormal_norm_clr[wlx_lognormal_norm_clr[,"p"] < 0.05, ], wlx_lognormal_norm_clr[wlx_lognormal_norm_clr[,"p"] < 0.05, ][1] != "NA")[1])
wlx_lognormal_norm_clr_outs
#css
wlx_lognormal_norm_css<- read.csv("/Users/jf/Desktop/BA/Analysis_SparseDOSSA/outputs/wilcoxon_output/wlx_pvalues_lognormal_norm_css.csv")
wlx_lognormal_norm_css_outs = as.list(subset(wlx_lognormal_norm_css[wlx_lognormal_norm_css[,"p"] < 0.05, ], wlx_lognormal_norm_css[wlx_lognormal_norm_css[,"p"] < 0.05, ][1] != "NA")[1])
wlx_lognormal_norm_css_outs
#tss
wlx_lognormal_norm_tss<- read.csv("/Users/jf/Desktop/BA/Analysis_SparseDOSSA/outputs/wilcoxon_output/wlx_pvalues_lognormal_norm_tss.csv")
wlx_lognormal_norm_tss_outs = as.list(subset(wlx_lognormal_norm_tss[wlx_lognormal_norm_tss[,"p"] < 0.05, ], wlx_lognormal_norm_tss[wlx_lognormal_norm_tss[,"p"] < 0.05, ][1] != "NA")[1])
wlx_lognormal_norm_tss_outs
#relab
wlx_lognormal_norm_relab<- read.csv("/Users/jf/Desktop/BA/Analysis_SparseDOSSA/outputs/wilcoxon_output/wlx_pvalues_lognormal_norm_relab.csv")
wlx_lognormal_norm_relab_outs = as.list(subset(wlx_lognormal_norm_relab[wlx_lognormal_norm_relab[,"p"] < 0.05, ], wlx_lognormal_norm_relab[wlx_lognormal_norm_relab[,"p"] < 0.05, ][1] != "NA")[1])
wlx_lognormal_norm_relab_outs
#spike
#clr
wlx_logspike_norm_clr<- read.csv("/Users/jf/Desktop/BA/Analysis_SparseDOSSA/outputs/wilcoxon_output/wlx_pvalues_logspike_norm_clr.csv")
wlx_logspike_norm_clr_outs = as.list(subset(wlx_logspike_norm_clr[wlx_logspike_norm_clr[,"p"] < 0.05, ], wlx_logspike_norm_clr[wlx_logspike_norm_clr[,"p"] < 0.05, ][1] != "NA")[1])
wlx_logspike_norm_clr_outs
#css
wlx_logspike_norm_css<- read.csv("/Users/jf/Desktop/BA/Analysis_SparseDOSSA/outputs/wilcoxon_output/wlx_pvalues_logspike_norm_css.csv")
wlx_logspike_norm_css_outs = as.list(subset(wlx_logspike_norm_css[wlx_logspike_norm_css[,"p"] < 0.05, ], wlx_logspike_norm_css[wlx_logspike_norm_css[,"p"] < 0.05, ][1] != "NA")[1])
wlx_logspike_norm_css_outs
#tss
wlx_logspike_norm_tss<- read.csv("/Users/jf/Desktop/BA/Analysis_SparseDOSSA/outputs/wilcoxon_output/wlx_pvalues_logspike_norm_tss.csv")
wlx_logspike_norm_tss_outs = as.list(subset(wlx_logspike_norm_tss[wlx_logspike_norm_tss[,"p"] < 0.05, ], wlx_logspike_norm_tss[wlx_logspike_norm_tss[,"p"] < 0.05, ][1] != "NA")[1])
wlx_logspike_norm_tss_outs
#relab
wlx_logspike_norm_relab<- read.csv("/Users/jf/Desktop/BA/Analysis_SparseDOSSA/outputs/wilcoxon_output/wlx_pvalues_logspike_norm_relab.csv")
wlx_logspike_norm_relab_outs = as.list(subset(wlx_logspike_norm_relab[wlx_logspike_norm_relab[,"p"] < 0.05, ], wlx_logspike_norm_relab[wlx_logspike_norm_relab[,"p"] < 0.05, ][1] != "NA")[1])
wlx_logspike_norm_relab_outs

#ttest
#clr
ttest_lognormal_norm_clr<- read.csv("/Users/jf/Desktop/BA/Analysis_SparseDOSSA/outputs/ttest_output/tt_pvalues_lognormal_norm_clr.csv")
ttest_lognormal_norm_clr_outs = as.list(subset(wlx_lognormal_norm_clr[wlx_lognormal_norm_clr[,"p"] < 0.05, ], wlx_lognormal_norm_clr[wlx_lognormal_norm_clr[,"p"] < 0.05, ][1] != "NA")[1])
ttest_lognormal_norm_clr_outs
#css
ttest_lognormal_norm_css<- read.csv("/Users/jf/Desktop/BA/Analysis_SparseDOSSA/outputs/ttest_output/tt_pvalues_lognormal_norm_css.csv")
ttest_lognormal_norm_css_outs = as.list(subset(wlx_lognormal_norm_css[wlx_lognormal_norm_css[,"p"] < 0.05, ], wlx_lognormal_norm_css[wlx_lognormal_norm_css[,"p"] < 0.05, ][1] != "NA")[1])
ttest_lognormal_norm_css_outs
#tss
ttest_lognormal_norm_tss<- read.csv("/Users/jf/Desktop/BA/Analysis_SparseDOSSA/outputs/ttest_output/tt_pvalues_lognormal_norm_tss.csv")
ttest_lognormal_norm_tss_outs = as.list(subset(wlx_lognormal_norm_tss[wlx_lognormal_norm_tss[,"p"] < 0.05, ], wlx_lognormal_norm_tss[wlx_lognormal_norm_tss[,"p"] < 0.05, ][1] != "NA")[1])
ttest_lognormal_norm_tss_outs
#relab
ttest_lognormal_norm_relab<- read.csv("/Users/jf/Desktop/BA/Analysis_SparseDOSSA/outputs/ttest_output/tt_pvalues_lognormal_norm_relab.csv")
ttest_lognormal_norm_relab_outs = as.list(subset(wlx_lognormal_norm_relab[wlx_lognormal_norm_relab[,"p"] < 0.05, ], wlx_lognormal_norm_relab[wlx_lognormal_norm_relab[,"p"] < 0.05, ][1] != "NA")[1])
ttest_lognormal_norm_relab_outs

#ttest spike
#clr
ttest_logspike_norm_clr<- read.csv("/Users/jf/Desktop/BA/Analysis_SparseDOSSA/outputs/ttest_output/tt_pvalues_logspike_norm_clr.csv")
ttest_logspike_norm_clr_outs = as.list(subset(wlx_logspike_norm_clr[wlx_logspike_norm_clr[,"p"] < 0.05, ], wlx_logspike_norm_clr[wlx_logspike_norm_clr[,"p"] < 0.05, ][1] != "NA")[1])
ttest_logspike_norm_clr_outs
#css
ttest_logspike_norm_css<- read.csv("/Users/jf/Desktop/BA/Analysis_SparseDOSSA/outputs/ttest_output/tt_pvalues_logspike_norm_css.csv")
ttest_logspike_norm_css_outs = as.list(subset(wlx_logspike_norm_css[wlx_logspike_norm_css[,"p"] < 0.05, ], wlx_logspike_norm_css[wlx_logspike_norm_css[,"p"] < 0.05, ][1] != "NA")[1])
ttest_logspike_norm_css_outs
#tss
ttest_logspike_norm_tss<- read.csv("/Users/jf/Desktop/BA/Analysis_SparseDOSSA/outputs/ttest_output/tt_pvalues_logspike_norm_tss.csv")
ttest_logspike_norm_tss_outs = as.list(subset(wlx_logspike_norm_tss[wlx_logspike_norm_tss[,"p"] < 0.05, ], wlx_logspike_norm_tss[wlx_logspike_norm_tss[,"p"] < 0.05, ][1] != "NA")[1])
ttest_logspike_norm_tss_outs
#relab
ttest_logspike_norm_relab<- read.csv("/Users/jf/Desktop/BA/Analysis_SparseDOSSA/outputs/ttest_output/tt_pvalues_logspike_norm_relab.csv")
ttest_logspike_norm_relab_outs = as.list(subset(wlx_logspike_norm_relab[wlx_logspike_norm_relab[,"p"] < 0.05, ], wlx_logspike_norm_relab[wlx_logspike_norm_relab[,"p"] < 0.05, ][1] != "NA")[1])
ttest_logspike_norm_relab_outs



#aldex
#lognorm welchs test
aldex2_lognormal_wel<- read.csv("/Users/jf/Desktop/BA/Analysis_SparseDOSSA/outputs/aldex2_output/aldex2_lognormal_welch.csv")
aldex2_lognormal_wel_outs = as.list(subset(aldex2_lognormal_wel[aldex2_lognormal_wel[,"we.ep"] < 0.05, ], aldex2_lognormal_wel[aldex2_lognormal_wel[,"we.ep"] < 0.05, ][1] != "NA")[2])
aldex2_lognormal_wel_outs
#lognorm wilcoxon
aldex2_lognormal_wil<- read.csv("/Users/jf/Desktop/BA/Analysis_SparseDOSSA/outputs/aldex2_output/aldex2_lognormal_wil.csv")
aldex2_lognormal_wil_outs = as.list(subset(aldex2_lognormal_wil[aldex2_lognormal_wil[,"wi.ep"] < 0.05, ], aldex2_lognormal_wil[aldex2_lognormal_wil[,"wi.ep"] < 0.05, ][1] != "NA")[2])
aldex2_lognormal_wil_outs
#spike welchs test
aldex2_spike_wel<- read.csv("/Users/jf/Desktop/BA/Analysis_SparseDOSSA/outputs/aldex2_output/aldex2_logspike_welch.csv")
aldex2_logspike_wel_outs = as.list(subset(aldex2_spike_wel[aldex2_spike_wel[,"we.ep"] < 0.05, ], aldex2_spike_wel[aldex2_spike_wel[,"we.ep"] < 0.05, ][1] != "NA")[2])
aldex2_logspike_wel_outs
#spike wilcoxon
aldex2_spike_wil<- read.csv("/Users/jf/Desktop/BA/Analysis_SparseDOSSA/outputs/aldex2_output/aldex2_logspike_wil.csv")
aldex2_logspike_wil_outs = as.list(subset(aldex2_spike_wil[aldex2_spike_wil[,"wi.ep"] < 0.05, ], aldex2_spike_wil[aldex2_spike_wil[,"wi.ep"] < 0.05, ][1] != "NA")[2])
aldex2_logspike_wil_outs


#lognormal
#aldex vs ttest
intersect_ttest_lognormal_norm_css_aldex_wil = intersect(ttest_lognormal_norm_css_outs,aldex2_lognormal_wil_outs)
intersect_ttest_lognormal_norm_css_aldex_wel = intersect(ttest_lognormal_norm_css_outs,aldex2_lognormal_wel_outs)
intersect_ttest_lognormal_norm_tss_aldex_wil = intersect(ttest_lognormal_norm_tss_outs,aldex2_lognormal_wil_outs)
intersect_ttest_lognormal_norm_tss_aldex_wel = intersect(ttest_lognormal_norm_tss_outs,aldex2_lognormal_wel_outs)
intersect_ttest_lognormal_norm_clr_aldex_wil = intersect(ttest_lognormal_norm_clr_outs,aldex2_lognormal_wil_outs)
intersect_ttest_lognormal_norm_clr_aldex_wel = intersect(ttest_lognormal_norm_clr_outs,aldex2_lognormal_wel_outs)
intersect_ttest_lognormal_norm_relab_aldex_wil = intersect(ttest_lognormal_norm_relab_outs,aldex2_lognormal_wil_outs)
intersect_ttest_lognormal_norm_relab_aldex_wel = intersect(ttest_lognormal_norm_relab_outs,aldex2_lognormal_wel_outs)

intersect_ttest_lognormal_norm_css_aldex_wil 
intersect_ttest_lognormal_norm_css_aldex_wel 
intersect_ttest_lognormal_norm_tss_aldex_wil 
intersect_ttest_lognormal_norm_tss_aldex_wel
intersect_ttest_lognormal_norm_clr_aldex_wil 
intersect_ttest_lognormal_norm_clr_aldex_wel 
intersect_ttest_lognormal_norm_relab_aldex_wil
intersect_ttest_lognormal_norm_relab_aldex_wel

#aldex vs wilcoxon
intersect_wlx_lognormal_norm_css_aldex_wil = intersect(wlx_lognormal_norm_css_outs,aldex2_lognormal_wil_outs)
intersect_wlx_lognormal_norm_css_aldex_wel = intersect(wlx_lognormal_norm_css_outs,aldex2_lognormal_wel_outs)
intersect_wlx_lognormal_norm_tss_aldex_wil = intersect(wlx_lognormal_norm_tss_outs,aldex2_lognormal_wil_outs)
intersect_wlx_lognormal_norm_tss_aldex_wel = intersect(wlx_lognormal_norm_tss_outs,aldex2_lognormal_wel_outs)
intersect_wlx_lognormal_norm_clr_aldex_wil = intersect(wlx_lognormal_norm_clr_outs,aldex2_lognormal_wil_outs)
intersect_wlx_lognormal_norm_clr_aldex_wel = intersect(wlx_lognormal_norm_clr_outs,aldex2_lognormal_wel_outs)
intersect_wlx_lognormal_norm_relab_aldex_wil = intersect(wlx_lognormal_norm_relab_outs,aldex2_lognormal_wil_outs)
intersect_wlx_lognormal_norm_relab_aldex_wel = intersect(wlx_lognormal_norm_relab_outs,aldex2_lognormal_wel_outs)

intersect_wlx_lognormal_norm_css_aldex_wil
intersect_wlx_lognormal_norm_css_aldex_wel
intersect_wlx_lognormal_norm_tss_aldex_wil
intersect_wlx_lognormal_norm_tss_aldex_wel
intersect_wlx_lognormal_norm_clr_aldex_wil
intersect_wlx_lognormal_norm_clr_aldex_wel
intersect_wlx_lognormal_norm_relab_aldex_wil
intersect_wlx_lognormal_norm_relab_aldex_wel

#ttest vs wilcoxon 
intersect_ttest_lognormal_norm_css_wlx_lognormal_norm_css = intersect(ttest_lognormal_norm_css_outs,wlx_lognormal_norm_css_outs)
intersect_ttest_lognormal_norm_css_wlx_lognormal_norm_tss= intersect(ttest_lognormal_norm_css_outs,wlx_lognormal_norm_tss_outs)
intersect_ttest_lognormal_norm_css_wlx_lognormal_norm_clr = intersect(ttest_lognormal_norm_css_outs,wlx_lognormal_norm_clr_outs)
intersect_ttest_lognormal_norm_css_wlx_lognormal_norm_relab = intersect(ttest_lognormal_norm_css_outs,wlx_lognormal_norm_relab_outs)
intersect_ttest_lognormal_norm_tss_wlx_lognormal_norm_css = intersect(ttest_lognormal_norm_tss_outs,wlx_lognormal_norm_css_outs)
intersect_ttest_lognormal_norm_tss_wlx_lognormal_norm_tss= intersect(ttest_lognormal_norm_tss_outs,wlx_lognormal_norm_tss_outs)
intersect_ttest_lognormal_norm_tss_wlx_lognormal_norm_clr = intersect(ttest_lognormal_norm_tss_outs,wlx_lognormal_norm_clr_outs)
intersect_ttest_lognormal_norm_tss_wlx_lognormal_norm_relab = intersect(ttest_lognormal_norm_tss_outs,wlx_lognormal_norm_relab_outs)
intersect_ttest_lognormal_norm_clr_wlx_lognormal_norm_css =  intersect(ttest_lognormal_norm_clr_outs,wlx_lognormal_norm_css_outs)
intersect_ttest_lognormal_norm_clr_wlx_lognormal_norm_tss= intersect(ttest_lognormal_norm_clr_outs,wlx_lognormal_norm_tss_outs)
intersect_ttest_lognormal_norm_clr_wlx_lognormal_norm_clr = intersect(ttest_lognormal_norm_clr_outs,wlx_lognormal_norm_clr_outs)
intersect_ttest_lognormal_norm_clr_wlx_lognormal_norm_relab = intersect(ttest_lognormal_norm_clr_outs,wlx_lognormal_norm_relab_outs)
intersect_ttest_lognormal_norm_relab_wlx_lognormal_norm_css = intersect(ttest_lognormal_norm_relab_outs,wlx_lognormal_norm_css_outs)
intersect_ttest_lognormal_norm_relab_wlx_lognormal_norm_tss= intersect(ttest_lognormal_norm_relab_outs,wlx_lognormal_norm_tss_outs)
intersect_ttest_lognormal_norm_relab_wlx_lognormal_norm_clr = intersect(ttest_lognormal_norm_relab_outs,wlx_lognormal_norm_clr_outs)
intersect_ttest_lognormal_norm_relab_wlx_lognormal_norm_relab = intersect(ttest_lognormal_norm_relab_outs,wlx_lognormal_norm_relab_outs)

intersect_ttest_lognormal_norm_css_wlx_lognormal_norm_css
intersect_ttest_lognormal_norm_css_wlx_lognormal_norm_tss
intersect_ttest_lognormal_norm_css_wlx_lognormal_norm_clr
intersect_ttest_lognormal_norm_css_wlx_lognormal_norm_relab
intersect_ttest_lognormal_norm_tss_wlx_lognormal_norm_css
intersect_ttest_lognormal_norm_tss_wlx_lognormal_norm_tss
intersect_ttest_lognormal_norm_tss_wlx_lognormal_norm_clr 
intersect_ttest_lognormal_norm_tss_wlx_lognormal_norm_relab
intersect_ttest_lognormal_norm_clr_wlx_lognormal_norm_css
intersect_ttest_lognormal_norm_clr_wlx_lognormal_norm_tss
intersect_ttest_lognormal_norm_clr_wlx_lognormal_norm_clr
intersect_ttest_lognormal_norm_clr_wlx_lognormal_norm_relab
intersect_ttest_lognormal_norm_relab_wlx_lognormal_norm_css
intersect_ttest_lognormal_norm_relab_wlx_lognormal_norm_tss
intersect_ttest_lognormal_norm_relab_wlx_lognormal_norm_clr
intersect_ttest_lognormal_norm_relab_wlx_lognormal_norm_relab

#logspike
#aldex vs ttest
intersect_ttest_logspike_norm_css_aldex_wil = intersect(ttest_logspike_norm_css_outs,aldex2_logspike_wil_outs)
intersect_ttest_logspike_norm_css_aldex_wel = intersect(ttest_logspike_norm_css_outs,aldex2_logspike_wel_outs)
intersect_ttest_logspike_norm_tss_aldex_wil = intersect(ttest_logspike_norm_tss_outs,aldex2_logspike_wil_outs)
intersect_ttest_logspike_norm_tss_aldex_wel = intersect(ttest_logspike_norm_tss_outs,aldex2_logspike_wel_outs)
intersect_ttest_logspike_norm_clr_aldex_wil = intersect(ttest_logspike_norm_clr_outs,aldex2_logspike_wil_outs)
intersect_ttest_logspike_norm_clr_aldex_wel = intersect(ttest_logspike_norm_clr_outs,aldex2_logspike_wel_outs)
intersect_ttest_logspike_norm_relab_aldex_wil = intersect(ttest_logspike_norm_relab_outs,aldex2_logspike_wil_outs)
intersect_ttest_logspike_norm_relab_aldex_wel = intersect(ttest_logspike_norm_relab_outs,aldex2_logspike_wel_outs)
intersect_ttest_logspike_norm_css_aldex_wil
intersect_ttest_logspike_norm_css_aldex_wel
intersect_ttest_logspike_norm_tss_aldex_wil
intersect_ttest_logspike_norm_tss_aldex_wel
intersect_ttest_logspike_norm_clr_aldex_wil
intersect_ttest_logspike_norm_clr_aldex_wel
intersect_ttest_logspike_norm_relab_aldex_wel
intersect_ttest_logspike_norm_relab_aldex_wil

#aldex vs wilcoxon
intersect_wlx_logspike_norm_css_aldex_wil = intersect(wlx_logspike_norm_css_outs,aldex2_logspike_wil_outs)
intersect_wlx_logspike_norm_css_aldex_wel = intersect(wlx_logspike_norm_css_outs,aldex2_logspike_wel_outs)
intersect_wlx_logspike_norm_tss_aldex_wil = intersect(wlx_logspike_norm_tss_outs,aldex2_logspike_wil_outs)
intersect_wlx_logspike_norm_tss_aldex_wel = intersect(wlx_logspike_norm_tss_outs,aldex2_logspike_wel_outs)
intersect_wlx_logspike_norm_clr_aldex_wil = intersect(wlx_logspike_norm_clr_outs,aldex2_logspike_wil_outs)
intersect_wlx_logspike_norm_clr_aldex_wel = intersect(wlx_logspike_norm_clr_outs,aldex2_logspike_wel_outs)
intersect_wlx_logspike_norm_relab_aldex_wil = intersect(wlx_logspike_norm_relab_outs,aldex2_logspike_wil_outs)
intersect_wlx_logspike_norm_relab_aldex_wel = intersect(wlx_logspike_norm_relab_outs,aldex2_logspike_wel_outs)
intersect_wlx_logspike_norm_css_aldex_wil
intersect_wlx_logspike_norm_css_aldex_wel
intersect_wlx_logspike_norm_tss_aldex_wil
intersect_wlx_logspike_norm_tss_aldex_wel
intersect_wlx_logspike_norm_clr_aldex_wil
intersect_wlx_logspike_norm_clr_aldex_wel
intersect_wlx_logspike_norm_relab_aldex_wil 
intersect_wlx_logspike_norm_relab_aldex_wel


#ttest vs wilcoxon 
intersect_ttest_logspike_norm_css_wlx_logspike_norm_css = intersect(ttest_logspike_norm_css_outs,wlx_logspike_norm_css_outs)
intersect_ttest_logspike_norm_css_wlx_logspike_norm_tss= intersect(ttest_logspike_norm_css_outs,wlx_logspike_norm_tss_outs)
intersect_ttest_logspike_norm_css_wlx_logspike_norm_clr = intersect(ttest_logspike_norm_css_outs,wlx_logspike_norm_clr_outs)
intersect_ttest_logspike_norm_css_wlx_logspike_norm_relab = intersect(ttest_logspike_norm_css_outs,wlx_logspike_norm_relab_outs)
intersect_ttest_logspike_norm_tss_wlx_logspike_norm_css = intersect(ttest_logspike_norm_tss_outs,wlx_logspike_norm_css_outs)
intersect_ttest_logspike_norm_tss_wlx_logspike_norm_tss= intersect(ttest_logspike_norm_tss_outs,wlx_logspike_norm_tss_outs)
intersect_ttest_logspike_norm_tss_wlx_logspike_norm_clr = intersect(ttest_logspike_norm_tss_outs,wlx_logspike_norm_clr_outs)
intersect_ttest_logspike_norm_tss_wlx_logspike_norm_relab = intersect(ttest_logspike_norm_tss_outs,wlx_logspike_norm_relab_outs)
intersect_ttest_logspike_norm_clr_wlx_logspike_norm_css =  intersect(ttest_logspike_norm_clr_outs,wlx_logspike_norm_css_outs)
intersect_ttest_logspike_norm_clr_wlx_logspike_norm_tss= intersect(ttest_logspike_norm_clr_outs,wlx_logspike_norm_tss_outs)
intersect_ttest_logspike_norm_clr_wlx_logspike_norm_clr = intersect(ttest_logspike_norm_clr_outs,wlx_logspike_norm_clr_outs)
intersect_ttest_logspike_norm_clr_wlx_logspike_norm_relab = intersect(ttest_logspike_norm_clr_outs,wlx_logspike_norm_relab_outs)
intersect_ttest_logspike_norm_relab_wlx_logspike_norm_css = intersect(ttest_logspike_norm_relab_outs,wlx_logspike_norm_css_outs)
intersect_ttest_logspike_norm_relab_wlx_logspike_norm_tss= intersect(ttest_logspike_norm_relab_outs,wlx_logspike_norm_tss_outs)
intersect_ttest_logspike_norm_relab_wlx_logspike_norm_clr = intersect(ttest_logspike_norm_relab_outs,wlx_logspike_norm_clr_outs)
intersect_ttest_logspike_norm_relab_wlx_logspike_norm_relab = intersect(ttest_logspike_norm_relab_outs,wlx_logspike_norm_relab_outs)
intersect_ttest_logspike_norm_css_wlx_logspike_norm_css
intersect_ttest_logspike_norm_css_wlx_logspike_norm_tss
intersect_ttest_logspike_norm_css_wlx_logspike_norm_clr
intersect_ttest_logspike_norm_css_wlx_logspike_norm_relab
intersect_ttest_logspike_norm_tss_wlx_logspike_norm_css
intersect_ttest_logspike_norm_tss_wlx_logspike_norm_tss
intersect_ttest_logspike_norm_tss_wlx_logspike_norm_clr
intersect_ttest_logspike_norm_tss_wlx_logspike_norm_relab
intersect_ttest_logspike_norm_clr_wlx_logspike_norm_css
intersect_ttest_logspike_norm_clr_wlx_logspike_norm_tss
intersect_ttest_logspike_norm_clr_wlx_logspike_norm_clr
intersect_ttest_logspike_norm_clr_wlx_logspike_norm_relab
intersect_ttest_logspike_norm_relab_wlx_logspike_norm_css
intersect_ttest_logspike_norm_relab_wlx_logspike_norm_tss
intersect_ttest_logspike_norm_relab_wlx_logspike_norm_clr 
intersect_ttest_logspike_norm_relab_wlx_logspike_norm_relab 





