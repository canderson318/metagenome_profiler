
#import "utils.typ": *
#import "@preview/zebraw:0.6.1": *
#show: zebraw

////////////////////
// Document Settings
////////////////////
#let fntsz = 11pt
#let after = fntsz * 1
#set page(margin: .5in, width: 8.5in, height: 11in, numbering: "1")
#set columns(gutter: fntsz)
#set text(font: "Arial", size: fntsz)
#set par(leading: fntsz, spacing: after,justify: true)

////////////////
// Set Defaults
////////////////

#show math.equation: set text(size: 12pt)

#show title: it => {
  set align(center + horizon)
  set text(size: 20pt, weight: "bold")
  v(-5em)
  it
}

#show outline: it => {
  show heading: set text(size: 15pt)
  show outline.entry: it => pad(left: 1em, right: 0em, it)
  it
}

#show link: it =>{
    set text(blue)
    it
}

#show ref: it =>{
    if it.element != none{
        set text(baseline: 0em, size: 1em, weight:"bold")
        it
    }else{
        set text(baseline: -.11em, size: .8em)
        it
    }
}

#show heading: it=>{
    block(below:after)[#it]
}

#show figure.where(kind: "code"): it => {
    // let col = rgb("#2746e397")
    // block(fill: none, stroke: col, inset: 8pt,below:5pt, radius: 20pt, width: 100%)[
    //     #it.body
    // ]
    block[
        #it.body
    ]
    v(-1em)
    it.caption
}

#let fig_width = 200pt
#show figure.caption: it =>{
    let sz = fntsz - 3pt
    set text(sz)
    set par(leading:sz/2)
    set align(left)
    h(10pt)
    box(width: fig_width )[
          *#it.supplement #it.counter.display(it.numbering)#it.separator* #it.body
    ]
}


#set heading(numbering:none)

#let theme_col= rgb("#8c2fc185")

#show "XXX": text.with(fill:red)

#show "GH": box[
    #box(
        link("https://github.com/canderson318/research-methodology-programming-assignment", image("../github.svg", height: .8em))
    )
]

#show heading.where(level:1): it =>{
    set text(fntsz+1pt,weight:"bold")
    it
}

#show heading.where(level:2): it =>{
    set text(fntsz,weight:"bold")
    it
}

#let fmt-num(n, decimals: 2, format: "f") = {
    if format == "e" {
      let exp = if n == 0 { 0 } else { calc.floor(calc.log(calc.abs(n), base: 10)) }
      let mantissa = n / calc.pow(10, exp)
      let m = str(calc.round(mantissa, digits: decimals))
      let parts = m.split(".")
      let frac = if parts.len() > 1 { parts.at(1) } else { "" }
      while frac.len() < decimals { frac = frac + "0" }
      parts.at(0) + "." + frac + "e+" + str(exp)
    } else {
      let s = str(calc.round(n, digits: decimals))
      let parts = s.split(".")
      let whole = parts.at(0)
      let frac = if parts.len() > 1 { parts.at(1) } else { "" }
      while frac.len() < decimals { frac = frac + "0" }
      let chars = whole.clusters().rev()
      let grouped = ()
      for (i, c) in chars.enumerate() {
        if i > 0 and calc.rem(i, 3) == 0 { grouped.push(",") }
        grouped.push(c)
      }
      grouped.rev().join("") + if decimals > 0 { "." + frac } else { "" }
    }
  }

//////////////////////
//////////////////////
//////// BODY ////////
//////////////////////
//////////////////////

// Written report (25%): This maximum 4-page (11pt arial font, ½ inch margins) technical document with an introduction, methods, and results sections.
// The introduction gives any necessary background and describes and justifies your strategy.
// The method section defines the input, presents the algorithm, and details the expected output.
// The results contain experiments, results, and discussion.
// You should cite the appropriate scientific literature where appropriate.
// The list of citations does not count against the page limit.

#counter(page).update(0)
#set page(numbering:none)
#title[Metatranscriptomic Virome Characterization of _Rhinolophus malayanus_ Bat Feces]
#align(center)[By: Christian Anderson\ #today.display("[day] [month repr:long] [year]")]

#align(center,
    stack([#v(10em)#image("../rhinopholus.jpg", height: 300pt)],[#v(.5em)#h(20em)@rhinolophus])
)

// #align(bottom, outline())
#pagebreak()

#set page(numbering:"1")

#set page(columns:2)

= Introduction
// - covid pandemic
// - importance of environmental monitoring
// - introduction of metagenome sampling
//     - shotgun sequencing
//     - host + microbiome
// - example with bats
//     - cite @Zhou2020

Environmental monitoring is a critical component of pandemic readiness.
In December 2019 cases of viral pneumonia began to arise due to some unidentified microbial agent.

Not soon after it was identified to be a novel strain of coronavirus (CoV), SARSCoV-2, whose primary host is the horseshoe bat (_Rhinolophidae_) native to the Yunnan Provence of China where the first cases were reported @Han2023.

A year later, Researchers were able to utilize shotgun sequencing and metagenomic analysis on a pooled sample across 78 feces samples to identify yet another CoV strain (RmYN02) present in the feces of the _Rhinolophodus malayanus_ horseshoe bat @Zhou2020.

Their findings revealed how CoV evolves and adapts to new hosts.
Outbreaks will continue to occur where humans encroach on animal spaces, more diseases will continue to jump hosts from animals to humans.
By continuing to monitor the source of viral pathogens, we can better predict its future spread and infection rates.

= Methods
// - profile microbiome by removing host signal
// - look at community abundances
//     - what communities are present
//     - which dominate
//     - (why do abundances matter?)
// - need simple method for characterizing virome
// - example with bats
// - host alignment and read filtering
//     - caveat: run on subset reference index so results are examples
//     - expected runtime with full genome
// - phred filtering
// - viral database querying 
// - abundance analysis
// - read edit distance analysis

I developed a general method as a python package that uses host-reference read filtering, quality filtering, and viral read classification to enable researchers to characterize the viral communities present in a short-read metagenomic sample.
This package takes as input a fastq file containing the short-read sequenced reads and outputs a suite of curated analyses after quality filtering for host and low quality reads.


I first indexed the host genome using K = 31 kmers.
I then generated the same size kmers for each of the sample reads and compared each set of read kmers to the set of host kmers.
Reads where greater than 50% of kmers matched a host kmer I deemed host originated and removed.
Next I calculated the average per-read Phred score and excluded reads with a average Phred less than 20.

The remaining reads were queried against a local copy of the Kraken2 viral database (#link("https://benlangmead.github.io/aws-indexes/k2#kraken2--bracken-16s-rna")[source]) to classify each read taxonomically @Wood2019.
Kraken2 uses exact k-mer matching to assign each read to the lowest common ancestor (LCA) in the database phylogeny, reporting a rank (species, genus, family, etc.) reflecting the resolution of the match.
To evaluate the overall classification resolution across reads, I calculated the ratio of species low to high rank classification ratio: 

#align(center)[$"count"("Species classifcations") / "count"("Genus,Family,Order, etc.
classifications")$.]
\
A ratio less than 1 means there were more high-rank classifications than low, indicating the reads are not diverse on the species level.

To calculate abundances I simply divided the count of each species by the total number of classified species, and the same for family and genus classifications.

As a complementary analysis, I calculated the string distance between each read and every other read.
String distance here is defined as the minimum number of edits to change one string into another, aka Levenshtein distance (LV).
I computed the LV between reads of only a stratified sample of 10,000 reads because $#fmt-num(39087705, decimals: 0)^2$ comparisons would take too long.
I then used Multi-dimensional Scaling (MDS) on the read edit-distance matrix to reduced the distances down to just two dimensions for plotting.
On the same matrix, I additionally performed a PERMANOVA test.
The PERMANOVA test tests whether the distance structure can be explained by a grouping variable, in this case species.
Specifically, it tests the null hypothesis centroids and dispersions are the same across groups @Bakker2026.

= Results

Filtering out host reads reduced the count down from #fmt-num(39257492, decimals: 0) to #fmt-num(39110352, decimals: 0), and phred filtering a further #fmt-num(22647, decimals: 0) left #fmt-num(39087705, decimals: 0) 'clean' reads which would ideally include only microbiome reads.
However, because ```python filter_reads()``` would take an estimated $2.77 times 10^12$ hours against the full genome (@runtime), the host index was built against only the $2^"nd"$ and $3^"rd"$ chromosome scaffolds, so a small fraction of remaining reads may still be host-derived.

#figure(caption:[Projected runtime with increasing scaffold count.
Three actual times and a log-linear model with an SE band in blue shows a exponential relationship.
Purple point shows expected time to index and align against the full 28 chromosome scaffolds.],
    align(left, [#image("../../out/figs/runtime_versus_index_scaffold_number.png", width: fig_width)]), 
) <runtime>

Projected runtimes for each scaffold count were estimated from a log-linear model of three iterations at 1, 2, and 3 scaffolds.
To reference against the full 28 chromosome scaffolds it would take #fmt-num(2.77e12, decimals: 0) hours.
To save myself #fmt-num(calc.round((2.77e12 / (100 * 360 * 24))),decimals:0) lifetimes worth of time I had to reference the reads against only the $2^"nd"$ and $3^"rd"$ chromosome scaffolds.

// total reads         =   39,257,492  (wc -l out/phreds.txt)
// host removed        =      147,140  (wc -l out/reads_from_host.inds)
// after host filter   =   39,110,352  = 39257492 - 147140
// phred removed       =       22,647  (wc -l out/phred_hits.txt)
// naive clean         =   39,087,705  = 39110352 - 22647
// actual clean        =   39,087,705  (grep -c ">" out/filtered.fasta)
// overlap (host+phred)=            0  reads flagged by both filters
// host %              = 0.3748%
// phred %             = 0.0577%
// total removed %     = 0.4325%
// 
#figure(align(center, image("../../out/figs/phred_distr.png", width: fig_width)), caption: [Read Phred score distribution.
The density peaks sharply around a Phred score of around 36–37, with very little density below 20.]) <phred>

The Phred distribution (@phred) confirms high overall sequencing quality: only #fmt-num(22647, decimals: 0) reads (0.06%) fell below the threshold of 20.

286,284 of 39,257,492 reads were classifiable against the Kraken2 viral database.
The species-to-genus/family rank ratio was 38:1, with species-level (S-rank) classifications outnumbering higher-rank assignments.
Again, a ratio greater than 1 indicates that the reads are sufficiently diverse and distinct for Kraken2's LCA to resolve them down to the species level rather than collapsing them into higher levels.
Genus, Family, and Species abundances show _Studiervirinae_, _Sarbercovirus_, and _Staphylococcus phage Andhra_ as having the highest abundances, respectively (@family_ab,@genus_ab,@species_ab).
Note, the reads not classifiable to Species are classified based on the next best alignment (Genus, then Family, etc.), so these abundances should be interpreted independently.

#grid(columns: (1fr,1fr),
    [
        #figure(caption:[\ Family classification \ log-abundances.],
            align(left, 
                image("../../out/figs/family_abundance_bar.png", width: fig_width - 50pt)
            )
        )<family_ab>
    ],
    [
        #figure(caption:[\ Genus classification \ log-abundances.],

            align(right, 
                image("../../out/figs/genus_abundance_bar.png", width: fig_width - 50pt)
            )
        )<genus_ab>
    ]
)

#figure(caption:[Top 90% Species classification log-abundances.],
    align(center, 
        image("../../out/figs/species_abundance_bar.png", width: fig_width+20pt)
    )
)<species_ab>

// Staphylococcus phage Andhra
// Severe acute respiratory syndrome coronavirus 2
// Escherichia phage 500465-1
// Proteus phage VB_PmiS-Isfahan
// SARS coronavirus Tor2
// SSRNA phage SRR5466369_2
// Vibrio phage ValB1MD-2
// Enterobacteria phage P7
// Escherichia phage VB_EcoM-613R3
// Shamonda virus
// Escherichia phage DE3
// SSRNA phage SRR5466337_3
// Rhizobium phage RHph_11
// Klebsiella phage ST147-VIM1phi7.1
// Corynebacterium phage Darwin
// Escherichia phage 500465-2
// Escherichia phage Bp7
// Choristoneura fumiferana granulovirus
// Klebsiella phage 4LV2017
// Escherichia phage 4A7
// Escherichia phage 2H10
// Physalis rugose mosaic virus
// Stx2-converting phage Stx2a_WGPS2
// Bat hepatitis E virus
// Stx2-converting phage 1717
// Choristoneura rosaceana entomopoxvirus 'L'
// Escherichia phage 1H12

The virome is dominated by bacteriophages (Staphylococcus phage Andhra, multiple _Escherichia_ phages, Proteus phage VB_PmiS-Isfahan, Vibrio phage ValB1MD-2), which together make up the majority of classified reads when considering the 90th percentile of abundances.
Corona virus also stands out as the second and fifth most abundant species.

The read string-distance MDS plot reveals clear clustering by species, indicating that reads from the same viral species are more similar to one another than to reads from other species (@mds).
PERMANOVA on the pairwise LV distance matrix confirmed this structure (pseudo-F = 133.5, p = .001), rejecting the null hypothesis that centroids and dispersions are equal across groups.
Together these results indicate that viral species identity is the primary driver of read-sequence variation in the sample.

#figure(caption:[MDS reduction of read Levenshtein distances colored by classified Species.
],
    align(center, 
        image("../../out/figs/top_perc_read_str_dist_mds.png", width: fig_width+80pt)
    )
)<mds>


= Discussion

Horseshoe bats (_Rhinolophidae_) have been well established as coronavirus resivoirs @Han2023, and this study's abundance results reinforce that finding.

SARS-CoV-2 and SARS coronavirus Tor2 rank among the top five most abundant species detected (@species_ab), and when taken together, represent the largest fraction of the classified virome, replicating Zhou et al.'s @Zhou2020 main result where _Rhinopholus malayanus_ fecal samples carry pandemic related CoV strains.
Shotgun metagenomics proves a valid and useful tool for monitoring pathogen zoonosis.

and thereby validates the utility of shotgun metagenomics as tool for monitoring pathogen zoonosis.
// This is expected: fecal samples capture the gut bacterial community alongside its associated phages, and gut-tropic bacteriophages routinely dominate metagenomic surveys of fecal viromes @Zhou2020.
// Bat hepatitis E virus also appears, consistent with the broad hepatotropic viral diversity documented in _Rhinolophus_ species.
// The Shannon diversity of 1.69 across all species indicates a moderately diverse community, not dominated to the point of excluding minority taxa.

The read string-distance analysis provides orthogonal confirmation of viral classifications.
The MDS embedding (@mds) shows that reads assigned to the same species have strong read similarity in LV space, while reads from different species are more distant.
The PERMANOVA result (pseudo-F = 133.5, p = .001) indicates that species membership is a strong and statistically significant predictor of read-sequence similarity.
This provides more confidence in the abundance estimates as it indicates that viral assignments are not arbitrary and reads from the same species are similar.

= AI use statement
AI was used in following cases during code development.
- Helped me understand file IO byte offsets and how to find specific line's byte offsets using ```python open(file).tell()```.
- Helped me understand how the line pointer shifts each time ```python open(file).readline()``` is called.
- Helped me understand how to use an initializer with ```python Multiprocessing.Pool()``` to share an object globally across workers.
- Helped me understand why github actions testing was failing when pytests were working locally.



#pagebreak()
#bibliography("../bibliography.bib", title: "References")