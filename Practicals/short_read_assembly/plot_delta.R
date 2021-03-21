args <- commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (delta file)\n", call.=FALSE)
}

delta <- args[1]

library(tidyverse)
library(magrittr)
library(GenomicRanges)
library(scales)

readDelta <- function(deltafile){
  lines = scan(deltafile, 'a', sep='\n', quiet=TRUE)
  lines = lines[-1]
  lines.l = strsplit(lines, ' ')
  lines.len = lapply(lines.l, length) %>% as.numeric
  lines.l = lines.l[lines.len != 1]
  lines.len = lines.len[lines.len != 1]
  head.pos = which(lines.len == 4)
  head.id = rep(head.pos, c(head.pos[-1], length(lines.l)+1)-head.pos)
  mat = matrix(as.numeric(unlist(lines.l[lines.len==7])), 7)
  res = as.data.frame(t(mat[1:5,]))
  colnames(res) = c('rs','re','qs','qe','error')
  res$qid = unlist(lapply(lines.l[head.id[lines.len==7]], '[', 2))
  res$rid = unlist(lapply(lines.l[head.id[lines.len==7]], '[', 1)) %>% gsub('^>', '', .)
  res$strand = ifelse(res$qe-res$qs > 0, '+', '-')
  res
}
filterMum <- function(df, minlength=1000, flanks=1e4){
    coord = df %>% filter(abs(re-rs)>minlength) %>% group_by(qid, rid) %>%
        summarize(qsL=min(qs)-flanks, qeL=max(qe)+flanks, rs=median(rs)) %>%
        ungroup %>% arrange(desc(rs)) %>%
        mutate(qid=factor(qid, levels=unique(qid))) %>% select(-rs)
    merge(df, coord) %>% filter(qs>qsL, qe<qeL) %>%
        mutate(qid=factor(qid, levels=levels(coord$qid))) %>% select(-qsL, -qeL)
}
diagMum <- function(df){
    ## Find best qid order
    rid.o = df %>% group_by(qid, rid) %>% summarize(base=sum(abs(qe-qs)),
                                                    rs=weighted.mean(rs, abs(qe-qs))) %>%
        ungroup %>% arrange(desc(base)) %>% group_by(qid) %>% do(head(., 1)) %>%
        ungroup %>% arrange(desc(rid), desc(rs)) %>%
        mutate(qid=factor(qid, levels=unique(qid)))
    ## Find best qid strand
    major.strand = df %>% group_by(qid) %>%
        summarize(major.strand=ifelse(sum(sign(qe-qs)*abs(qe-qs))>0, '+', '-'),
                  maxQ=max(c(qe, qs)))
    merge(df, major.strand) %>% mutate(qs=ifelse(major.strand=='-', maxQ-qs, qs),
                                       qe=ifelse(major.strand=='-', maxQ-qe, qe),
                                       qid=factor(qid, levels=levels(rid.o$qid)))
}
maxSimilarityDisjoin <- function(df){
  ref.ir = GRanges('X', IRanges(df$rs, df$re), similarity=df$similarity)
  ## Efficient clean up of low similarity within high similarity
  step = 1
  while(step>0){
    largealign = ref.ir[head(order(rank(-ref.ir$similarity), rank(-width(ref.ir))),step*1000)]
    ol = findOverlaps(ref.ir, largealign, type='within') %>% as.data.frame %>%
        mutate(simW=ref.ir$similarity[queryHits],
               simL=largealign$similarity[subjectHits]) %>% filter(simW<simL)
    if(length(largealign) == length(ref.ir)){
      step = 0
    } else {
      step = step + 1
    }
    ref.ir = ref.ir[-ol$queryHits]
  }
  ## Disjoin and annotate with the max similarity
  ref.dj = disjoin(c(ref.ir, GRanges('X', IRanges(min(df$rs), max(df$rs)), similarity=0)))
  ol = findOverlaps(ref.ir, ref.dj) %>% as.data.frame %>%
      mutate(similarity=ref.ir$similarity[queryHits]) %>%
      group_by(subjectHits) %>% summarize(similarity=max(similarity))
  ref.dj$similarity = 0
  ref.dj$similarity[ol$subjectHits] = ol$similarity
  as.data.frame(ref.dj)
}

# Load the delta file
mumgp <- readDelta(delta)

# Filter for alignments >=1kbp
mumgp.filt <- filterMum(mumgp, minlength=1000)

# Find query sequence ordering
mumgp.filt.diag <- diagMum(mumgp.filt)

mumgp %<>% mutate(similarity=1-error/abs(qe-qs))
mumgp.filt %<>% mutate(similarity=1-error/abs(qe-qs))


pdf(file=paste(delta, ".pdf", sep=""), width=7, height=7)

# Plot reference coverage and contig identity plots
ggplot(mumgp, aes(x=rs, xend=re, y=similarity*100, yend=similarity*100)) + geom_segment() +
    theme_bw() + xlab('Reference Sequence [bp]') + ylab('Identity [%]') + ggtitle('All contigs') +
    ylim(50,100)
ggplot(mumgp.filt, aes(x=rs, xend=re, y=similarity*100, yend=similarity*100)) + geom_segment() +
    theme_bw() + xlab('Reference Sequence [bp]') + ylab('Identity [%]') + ggtitle('At least 1 kbp aligned') +
    ylim(50,100)

# Plot alignments
ggplot(mumgp.filt.diag, aes(x=rs, xend=re, y=qs, yend=qe, colour=strand)) +
    geom_segment() +
    geom_point(
      alpha = 0.5
    ) +
    theme_bw() + 
    facet_grid(
      qid~rid,
      scales = 'free',
      space  = 'free',
      switch = 'y'
    ) +
    theme(
      strip.text.y         = element_text(angle=180, size=5),
      strip.background     = element_blank(),
      legend.position      = c(.99,.01),
      legend.justification = c(1,0),
      axis.text.y          = element_blank(),
      axis.ticks.y         = element_blank()
    ) +
    scale_x_continuous(labels = comma) +
    xlab('Reference Sequence [bp]') +
    ylab('Assembly') +
    scale_colour_brewer(palette='Set1')


#mumgp.sim = maxSimilarityDisjoin(mumgp)
#mumgp.filt.sim <- maxSimilarityDisjoin(mumgp.filt)
#mumgp.filt.m <- rbind(mumgp.sim %>% mutate(filter='before'),
#		                          mumgp.filt.sim %>% mutate(filter='after'))
#
#mumgp.filt.m %>% filter(similarity==0) %>%
#    ggplot(aes(x=start, xend=end, y=filter, yend=filter)) +
#    geom_segment(size=10) +
#    theme_bw() +
#    xlab('Reference Sequence [bp]') +
#    ylab('filter') +
#    scale_colour_brewer(palette='Set1') +
#    ggtitle('Reference regions not covered')

dev.off()
