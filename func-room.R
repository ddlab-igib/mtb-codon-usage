read_fasta_and_compute_codon<-function(fa_file){
  ffn_file_content<-readLines(fa_file,warn = FALSE);
  ##
  idx_header<-grep('^>',ffn_file_content);
  all_fasta<-list();
  for(i in 1:length(idx_header)){
    if(i != length(idx_header)){
      cur_idx<-idx_header[i];
      nxt_idx<-idx_header[i+1];
      fasta_header<-ffn_file_content[cur_idx];
      fasta_seq<-paste(ffn_file_content[(cur_idx+1):(nxt_idx-1)],collapse = '');
      all_fasta[[i]]<-list(fasta_header=fasta_header,
                           fasta_seq=fasta_seq);
    }else{
      cur_idx<-idx_header[i];
      last_idx<-length(ffn_file_content);
      fasta_header<-ffn_file_content[cur_idx];
      fasta_seq<-paste(ffn_file_content[(cur_idx+1):last_idx],collapse = '');
      all_fasta[[i]]<-list(fasta_header=fasta_header,
                           fasta_seq=fasta_seq);
    }
  }
  my_criteria<-function(length){
    return(paste(".{",length,"}",sep=''));
  }
  get_codon<-function(s){
    output<-regmatches(s, gregexpr(my_criteria(length = 3), s))[[1]];
    ### Above function doesnot include any string of undefined length
    return(output);
  }
  ###
  get_freq<-function(input){
    output<-table(get_codon(s=input));
    return(output);
  }
  library(foreach);
  codon_freq_org<-foreach(i = 1:length(all_fasta), .combine = rbind)%do%{
    df_freq<-data.frame(
      get_freq(input=all_fasta[[i]]$fasta_seq),
      fastaHeader = all_fasta[[i]]$fasta_header,
      stringsAsFactors = FALSE);
    colnames(df_freq)[1]<-'codons';
    df_freq
  }
  #codon_freq_org$codons <- as.character(codon_freq_org$codons)
  codon_freq_org$codons <- toupper(codon_freq_org$codons)
  return(codon_freq_org)
}
#####
# extract from gff
getLocus <- function(x){
  gsub('>','',strsplit(x,split='\\|')[[1]][1])
}
getGene <- function(x){
  strsplit(x,split='\\|')[[1]][2]
}
getType <- function(x){
  strsplit(x,split='\\|')[[1]][3]
}
###########
# for summary
drawDonut <- function(x, y, prop,
                      stroke.col = c('#ffcc00ff','#2a7fffff'),
                      r = 1.5,
                      stroke.size = 0.5,
                      label){
  proT <- prop[1]
  getAng <- function(proT){
    mid <- (360*proT)/100
    angs <- list(c(mid,359+0.5)+90,
                 c(1,mid)+89
    )
  }
  angs <- getAng(proT)
  for(i in 1:length(angs)){
    opac.arc(x = x,
             y = y,
             r = r,
             fill = 'white',
             stroke.size = stroke.size,
             start.angle = angs[[i]][1],
             end.angle = angs[[i]][2],
             stroke.col = stroke.col[i]
    )
  }
  text(x = x,
       y = y,
       labels = label)
}
