# Scale pseudotime from 0 to 1
range01 <- function(x){(x-min(x))/(max(x)-min(x))}

umap_plot=function(s.object,plot_variable,pt.size,text_size,color_scheme=NULL,filename,dpi,transparency=NULL,type,legend_position="none",ver="v1"){
  if (is.character(color_scheme)) {
    p1=UMAPPlot(s.object,group.by=plot_variable,pt.size=pt.size)+theme_void(base_size = text_size)+theme(legend.position = legend_position)+scale_color_manual(values = color_scheme)
  } else {
    p1=UMAPPlot(s.object,group.by=plot_variable,pt.size=pt.size)+theme_void(base_size = text_size)+theme(legend.position = legend_position)
  }
  if (is.double(transparency)) {
    p1[[1]]$layers[[1]]$aes_params$alpha = transparency
  }
  if (type=="manuscript") {
    ggsave(filename = filename,path = paste("~/summary_time_series/figures",ver,"manuscript/tiffs/",sep = "/"), width = 6, height = 4, device='tiff', dpi=dpi,units = "in")
  } else {
    ggsave(filename = filename,path = paste("~/summary_time_series/figures",ver,"manuscript/tiffs/",sep = "/"), width = 6, height = 4, device='tiff', dpi=dpi,units = "in")
  }
}

feature_plot=function(s.object,gene,pt.size,text_size,filename,dpi,type,ver="v1"){
  FeaturePlot(s.object,gene,pt.size=pt.size)+scale_color_viridis(option = "B")+theme_void(base_size = text_size)+theme(legend.position = "right")
  if (type=="manuscript") {
    ggsave(filename = filename,path = paste("~/summary_time_series/figures",ver,"manuscript/tiffs/",sep = "/"), width = 7.5, height = 5.5, device='tiff', dpi=dpi,units = "in")
  } else {
    ggsave(filename = filename,path = paste("~/summary_time_series/figures",ver,"manuscript/tiffs/",sep = "/"), width = 7.5, height = 5.5, device='tiff', dpi=dpi,units = "in")
  }
}

pseudoplot=function(s.object,cell_order_variable,pt.size,text_size,filename,dpi,type){
  s.object$cell_order_variable_scaled=range01(s.object@meta.data[,cell_order_variable])
  FeaturePlot(s.object,"cell_order_variable_scaled",pt.size=pt.size)+scale_color_gradientn(colors = brewer.pal(11,'Spectral')[-6])+theme_void(base_size = text_size)
  if (type=="manuscript") {
    ggsave(filename = filename,path = "~/summary_time_series/figures/v1/manuscript/tiffs/", width = 7.5, height = 5.5, device='tiff', dpi=dpi,units = "in")
  } else {
    ggsave(filename = filename,path = "~/summary_time_series/figures/v1/supplement/tiffs/", width = 7.5, height = 5.5, device='tiff', dpi=dpi,units = "in")
  }
}

plot_hm_cs=function(s.object,ordered_gene_list,cell_order_variable,gene_split,show_row_names,col_fun=NULL){
  cell_order_slingshot=rownames(s.object@meta.data)[order(s.object@meta.data[,cell_order_variable],decreasing = F)]
  s.object=ScaleData(s.object,features = ordered_gene_list)
  plot_var=s.object@assays$RNA@scale.data[ordered_gene_list,cell_order_slingshot]
  ht=Heatmap(plot_var,show_column_names = FALSE,show_row_names = show_row_names,cluster_columns = FALSE, row_split = gene_split,cluster_row_slices = F,col=col_fun)
  ht=draw(ht)
  return(ht)
}

plot_gene_list=function(s.object,gene_list,cell_order_variable,text_size,type,filename,comment_num=NA){
  s.object=ScaleData(s.object,features = gene_list)
  x=t(s.object@assays$RNA@scale.data[gene_list,])
  x=data.frame(x)
  x$cells=rownames(x)
  x$pd=s.object@meta.data[,cell_order_variable]
  x$pd=range01(x$pd)
  x=gather(x,"gene","data",-c(cells,pd))
  x$gene=factor(x$gene,levels = gene_list)
  ggplot(x,aes(x=pd,y=data,color=gene))+geom_smooth(size=1.5)+theme_classic(base_size = text_size)+theme(legend.position = "none")+scale_color_manual(values = colorRampPalette(brewer.pal(n = 8,name = "Set1"))(length(gene_list)))+labs(x="Pseudotime",y="Scaled Expression")
  if (type=="manuscript") {
    ggsave(filename = filename,path = "~/summary_time_series/figures/v1/manuscript/tiffs/", width = 6.5, height = 4.5, device='tiff', dpi=dpi,units = "in")
  } else if (type=="supplement") {
    ggsave(filename = filename,path = "~/summary_time_series/figures/v1/supplement/tiffs/", width = 6.5, height = 4.5, device='tiff', dpi=dpi,units = "in")
  }
  else {
    ggsave(filename = filename,path = paste("~/summary_time_series/revisions/",comment_num,"/tiffs/",sep = ""), width = 6.5, height = 4.5, device='tiff', dpi=dpi,units = "in")
  }
}

plotgo=function(path_to_file,source,topterms,text_size,type,filename,dpi){
  if (source=="geneontology.org") {
    go_terms=read.delim(path_to_file, header=FALSE, comment.char="#")
    go_terms=go_terms[,c(1,3,ncol(go_terms))]
    colnames(go_terms)=c("Term","Entities","FDR")
    go_terms=subset(go_terms,subset=FDR<0.05)
    go_terms=subset(go_terms,subset=Entities<20)
    go_terms=go_terms[order(go_terms$FDR,decreasing = F),]
    ggplot(go_terms[1:topterms,],aes(x=reorder(Term, -log10(FDR)),y=-log10(FDR)))+geom_bar(stat = "identity",fill="white",color="navyblue",size=1.5)+coord_flip()+theme_classic(base_size=text_size)+xlab("GO Terms")
  }
  if (type=="manuscript") {
    ggsave(filename = filename,path = "~/summary_time_series/figures/v1/manuscript/tiffs/", width = 6, height = 4, device='tiff', dpi=dpi,units = "in")
  } else {
    ggsave(filename = filename,path = "~/summary_time_series/figures/v1/supplement/tiffs/", width = 6, height = 4, device='tiff', dpi=dpi,units = "in")
  }
}

barplotgenebymodule=function(gene_list,gene_cluster,lineage,text_size){
  if (lineage=="white"){
    plot.df=as.data.frame(table(gene_cluster[gene_list]))
    plot.df$module=NA
    plot.df$module[plot.df$Var1==0]="gradual_downregulation"
    plot.df$module[plot.df$Var1==1]="delayed_upregulation"
    plot.df$module[plot.df$Var1==2]="immediate_downregulation"
    plot.df$module[plot.df$Var1==3]="transient_upregulation"
    plot.df$module[plot.df$Var1==4]="gradual_upregulation"
    p=ggplot(plot.df,aes(x=reorder(module,Freq),y=Freq))+geom_bar(stat = "identity",fill="white",color="navyblue",size=1)+xlab("Module")+ylab("Count")+theme_classic(base_size = text_size)+coord_flip()
  } else{
    plot.df=as.data.frame(table(gene_cluster[gene_list]))
    plot.df$module=NA
    plot.df$module[plot.df$Var1==0]="delayed_upregulation"
    plot.df$module[plot.df$Var1==1]="immediate_downregulation"
    plot.df$module[plot.df$Var1==2]="transient_upregulation"
    plot.df$module[plot.df$Var1==3]="gradual_downregulation"
    plot.df$module[plot.df$Var1==4]="gradual_upregulation"
    p=ggplot(plot.df,aes(x=reorder(module,Freq),y=Freq))+geom_bar(stat = "identity",fill="white",color="navyblue",size=1)+xlab("Module")+ylab("Count")+theme_classic(base_size = text_size)+coord_flip()
  }
  return(p)
}

plotclusterbyday=function(s.object,text_size,color_scheme,filename,dpi,legend_position="top",ver="v1"){
  x=data.frame(round(table(s.object$batch.indices,s.object$cluster)/rowSums(table(s.object$batch.indices,s.object$cluster))*100,1))
  colnames(x)=c("Day","Cluster","Percent")
  ggplot(x,aes(x=Day,y=Percent,fill=Cluster))+geom_bar(stat = "identity",position = position_dodge(),color="black")+theme_classic(base_size = size)+scale_fill_manual(values = color_scheme)+theme(legend.position = legend_position)
  ggsave(filename = filename,path = paste("~/summary_time_series/figures",ver,"manuscript/tiffs/",sep = "/"), width = 6, height = 4, device='tiff', dpi=dpi,units = "in")
}
