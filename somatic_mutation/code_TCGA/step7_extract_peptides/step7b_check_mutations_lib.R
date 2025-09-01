
update_df_count <- function(df_counts, mut1_all, mut_data1, var_class){
  df_counts[Variant_Classification == var_class, 
            n_all := nrow(mut1_all)]
  
  df_counts[Variant_Classification == var_class, 
            n_mut5_no := nrow(mut_data1[group=="No"])]
  
  df_counts[Variant_Classification == var_class, 
            n_mut5_weak := nrow(mut_data1[group=="Weak"])]
  
  df_counts[Variant_Classification == var_class, 
            n_mut5_strong := nrow(mut_data1[group=="Strong"])]
  
  df_counts
}

get_DNA_ref_mut_plot <- function(DNA_ref_all, DNA_mut_all, 
                                 mut_data1, DNA_ref, DNA_mut, 
                                 freq_upper_limit = 1){
  tb0 = table(DNA_ref_all, DNA_mut_all)/length(DNA_ref_all)
  
  tb1 = list()
  for(g1 in unique(mut_data1$group)){
    w1 = which(mut_data1$group == g1)
    tb1[[g1]] = table(DNA_ref[w1], DNA_mut[w1])/length(DNA_ref[w1])
  }
  
  g1_list = list()
  
  g1_list[["All"]] = ggplot(as.data.frame(round(tb0, 2)), 
                            aes(x=DNA_ref_all, y=DNA_mut_all, fill=Freq)) +
    geom_tile() + ggtitle("All") + 
    geom_text(aes(label=Freq), color="black", size=3.5, na.rm = TRUE) + 
    scale_fill_gradient(low="white", high="red", limit=c(0,freq_upper_limit)) + 
    xlab("Reference allele") + ylab("Mutated allele") + 
    theme(legend.title = element_blank())
  
  for(grp1 in c("No", "Weak", "Strong")){
    g1_list[[grp1]] = 
      ggplot(as.data.frame(round(tb1[[grp1]], 2)), 
             aes(x=Var1, y=Var2, fill=Freq)) +
      geom_tile() + ggtitle(grp1) + 
      geom_text(aes(label=Freq), color="black", size=3.5, na.rm = TRUE) + 
      scale_fill_gradient(low="white", high="red", limit=c(0,freq_upper_limit)) + 
      xlab("Reference allele") + ylab("Mutated allele") + 
      theme(legend.title = element_blank())
  }
  
  g1_list
}


plot_AA_ref_mut_all <- function(AA_ref_all, AA_mut_all, legend_on_top = FALSE){
  
  df0 = as.data.frame(table(AA_ref_all, AA_mut_all))
  names(df0) = c("AA_ref", "AA_mut", "Freq")
  df0[["Prop"]] = df0$Freq/sum(df0$Freq)
  
  gg0 = ggplot(df0 %>% mutate(Prop = round(Prop, 2)), 
               aes(x=AA_ref, y=AA_mut, fill=Prop)) +
    geom_tile() + ggtitle("all mutations") + 
    geom_text(aes(label=Prop), color="black", size=3, na.rm = TRUE) + 
    scale_fill_gradient(low="white", high="red") + 
    xlab("Reference allele") + ylab("Mutated allele")
  
  if(legend_on_top){
    gg0 = gg0 + theme(legend.position = "top")
  }
  print(gg0)
}

plot_AA_ref_mut <- function(AA_ref_all, AA_mut_all, 
                            mut_data1, AA_ref, AA_mut, 
                            legend_on_top = FALSE){
  
  df0 = as.data.frame(table(AA_ref_all, AA_mut_all))
  names(df0) = c("AA_ref", "AA_mut", "Freq")
  df0[["Prop"]] = df0$Freq/sum(df0$Freq)
  
  df1 = list()
  for(g1 in c("No", "Weak", "Strong")){
    w1 = which(mut_data1$group == g1)
    tb_g1 = table(AA_ref[w1], AA_mut[w1])
    df_g1 = as.data.frame(tb_g1)
    names(df_g1) = c("AA_ref", "AA_mut", "Freq")
    df_g1[["Prop"]] = df_g1$Freq/sum(df_g1$Freq)
    df1[[g1]] = df_g1
    
    gg_g1 = ggplot(df_g1 %>% mutate(Prop = round(Prop, 2)), 
                   aes(x=AA_ref, y=AA_mut, fill=Prop)) +
      geom_tile() + ggtitle(g1) + 
      geom_text(aes(label=Prop), color="black", size=3, na.rm = TRUE) + 
      scale_fill_gradient(low="white", high="red") + 
      xlab("Reference allele") + ylab("Mutated allele")
    
    if(legend_on_top){
      gg_g1 = gg_g1 + theme(legend.position = "top")
    }
    print(gg_g1)
  }
  
  
  dfa = merge(df0, df1[["Strong"]], by=c("AA_ref", "AA_mut"), 
              suffixes = c("_All", "_Strong"))
  
  dfb = merge(df1[["Weak"]], df1[["No"]], by=c("AA_ref", "AA_mut"), 
              suffixes = c("_Weak", "_No"))
  
  dfa = merge(dfa, dfb, by=c("AA_ref", "AA_mut"))
  
  dim(dfa)
  dfa[1:2,]
  
  dfa = dfa %>% mutate(Ratio = Prop_Strong/Prop_All)
  dfa = dfa[order(dfa$Ratio, decreasing = TRUE),]
  
  dfa$odds_ratio = dfa$pvalue = rep(NA, nrow(dfa))

  for(i in 1:nrow(dfa)){
    x_i = matrix(NA, nrow=2, ncol=2)
    x_i[1,1] = sum(dfa$Freq_All[-i]) - sum(dfa$Freq_Strong[-i])
    x_i[1,2] = sum(dfa$Freq_Strong[-i])
    x_i[2,1] = dfa$Freq_All[i] - dfa$Freq_Strong[i]
    x_i[2,2] = dfa$Freq_Strong[i]
    fi = fisher.test(x_i)
    dfa$odds_ratio[i] = fi$estimate
    dfa$pvalue[i] = fi$p.value
  }
  
  gg_ratio = ggplot(dfa %>% mutate(Ratio = round(Ratio, 2)),  
                    aes(x=AA_ref, y=AA_mut, fill=-log10(pvalue))) +
    geom_tile() + ggtitle("potential neoantigen/all mutations") + 
    geom_text(aes(label=Ratio), color="black", size=3, na.rm = TRUE) + 
    scale_fill_gradient2(low="blue", high="red", mid="white", 
                         midpoint=2, na.value="white") + 
    xlab("Reference allele") + ylab("Mutated allele")
  
  if(legend_on_top){
    gg_ratio = gg_ratio + theme(legend.position = "top")
  }
  print(gg_ratio)
  
  dfa
}

plot_AA_prop_all_vs_strong <- function(df_sub){
  df_long = df_sub %>%
    pivot_longer(
      cols = c(Prop_All, Prop_Strong),
      names_to = "type",
      values_to = "proportion"
    )
  
  df_long = df_long %>%
    unite(mutation, AA_ref, AA_mut, sep = " > ")
  
  df_long = df_long %>%
    mutate(type = gsub("Prop_", "", type, fixed=TRUE))

  gg1 = ggplot(df_long, aes(x = mutation, y = proportion, fill = type)) + 
    geom_bar(stat = "identity", position = "dodge") +
    labs(x = "Mutation: Refeerence > Mutated", y = "Proportion", fill = "") +
    scale_fill_manual(values = c("All" = "skyblue", "Strong" = "orange"))
  
  print(gg1)
}

plot_AA_prop_all_vs_public <- function(df_sub){
  df_long = df_sub %>%
    pivot_longer(
      cols = c(Prop_All, Prop_Strong, Prop_Weak, Prop_No),
      names_to = "type",
      values_to = "proportion"
    )
  
  df_long = df_long %>%
    unite(mutation, AA_ref, AA_mut, sep = " > ")
  
  df_long = df_long %>%
    mutate(type = gsub("Prop_", "", type, fixed=TRUE))
  
  colors = c("skyblue", brewer.pal(n = 8, name = "Accent")[1:3])

  df_long$type = factor(df_long$type, levels = c("All", "No", "Weak", "Strong"))
  
  gg1 = ggplot(df_long, aes(x = mutation, y = proportion, fill = type)) + 
    geom_bar(stat = "identity", position = "dodge") +
    labs(x = "Mutation: Refeerence > Mutated", y = "Proportion", fill = "") +
    scale_fill_manual(values = colors)
  
  print(gg1)
}
