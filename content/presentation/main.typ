
#import "lib.typ": *

#let theme_col= rgb("5a0d87")
#let foot_col= red

#show footnote.entry: set text(foot_col,size:.8em)

#show footnote: it => (
    context {
      text(fill: foot_col, size: 1.2em)[#it]
    }
  )

// #show figure.caption: emph

#show figure.caption: it =>{
    let sz = .8em
    set text(sz)
    set par(leading:sz/2)
    set align(center)

    box(width: 500pt )[
          *#it.supplement #it.counter.display(it.numbering)#it.separator* #it.body
    ]
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

#show: typslides.with(
  ratio: "16-9",
  theme: theme_col,
  font: "Fira Sans",
  font-size: 20pt,
  link-style: "color",
  show-progress: true,
)

// ############ ############ ############ ############ ############
// ############ ############ ############ ############ ############
// ############ ############ ############ ############ ############
// ############ ############ ############ ############ ############

#front-slide(
  title: "Metagenome Profiler",
  authors: "Christian Anderson",
  info: today.display( "[day] [month repr:long] [year]"),
)

#slide(title:"Introduction")[
  #grid(columns: (1fr,.4fr), gutter: 20pt)[
    - SARS-CoV-2 emerged December 2019
      - primary host: horseshoe bat (_Rhinolophidae_), Yunnan, China @Han2023
    - Zhou et al. @Zhou2020 identified novel CoV strain (RmYN02) in _Rhinolophus malayanus_ feces via shotgun metagenomics 
    - Outbreaks will continue as humans encroach on animal spaces
    - Environmental monitoring is critical for pandemic readiness
    - Monitoring viral communities enables prediction of future spread
  ][
    #stack(dir:ttb,
      image("../rhinopholus.jpg", height: 200pt),
      { align(right)[@rhinolophus]}
    )
  ]
]

#slide(title:"Methods")[
  
  == Data

  *_Reads_*
  - Short-read shotgun metatranscriptomic sequencing of feces from _Rhinolophus malayanus_ bats @Zhou2020
  *_Reference Genome_*
  - _Rhinolophus ferrumequinum_ shotgun WGS assembly 
    - Closest one I could find to _R. malayanus_

  #pagebreak()

  // #grid(columns: (1fr, 1fr), gutter: 20pt)[
    #set align(horizon)

    ==  Pipeline
    + *Host read removal*: \ ~~~K=31 kmer index of host genome; reads with >50% kmer match removed
    + *Quality filtering*: \ ~~~reads with avg Phred < 20 removed
    + *Viral classification*: \ ~~~Kraken2 LCA exact k-mer matching against viral database @Wood2019
    + *Abundance analysis*: \ ~~~species/genus/family counts ÷ total classified
    + *Read edit distance analysis*: \ ~~~pairwise Levenshtein distance $arrow.r$ MDS + PERMANOVA
  // ]
  // [
  //   #set text(.8em)
  //   #set align(bottom)
  //   _*Classification resolution*_

  //   #text(.8em)[$ "count"("Species") / "count"("Genus, Family, ...") $]

  //   Ratio < 1 $arrow.r$ reads not diverse at species level

  //   _*String distance*_

  //   - Levenshtein distance on stratified sample of 10,000 reads
  //     - MDS $arrow.r$ 2D
  //     - PERMANOVA test: 
  //       - Does species membership explain distance structure?
  // ]
]


#slide(title:"Results: Filtering")[
  #set text(.9em)
  #grid(columns: (1fr, auto), gutter: 20pt)[
    - Started with #fmt-num(39257492, decimals:0) reads
    - full 28-chromosome run estimated at #fmt-num(2.77e12, decimals:0) hours (@runtime)
    #pad(left:40pt)[*$arrow.r.curve$* index built against only the $2^"nd"$ and $3^"rd"$ chromosome scaffolds]
    - Host filter (chr. 2 & 3) removed #fmt-num(147140, decimals:0) reads #strong[$arrow.r$ #fmt-num(39110352, decimals:0)]
  ][
    #figure(
      image("../../out/figs/runtime_versus_index_scaffold_number.png", width: 300pt),
      caption: [Projected runtime vs.\ scaffold count.]
    ) <runtime>
  ]
]

#slide(title:"Results: Read Quality")[
  #set align(top)
  #figure(
    align(center, image("../../out/figs/phred_distr.png", width: 400pt)),
    caption: [Read Phred score distribution. \ Density peaks around 36–37; very little density below 20.]
  ) <phred>

  #v(-.5em)
  - Only *#fmt-num(22647, decimals:0) reads (0.06%) fell below Phred threshold of 20* 
    - #text(theme_col, weight:"bold")[high overall sequencing quality.]
]

#slide(title:"Results: Abundances")[
  - *286,284* of #fmt-num(39257492, decimals:0) reads classified against Kraken2 viral database.
  - \#Species *:* \#Genus/Family ratio of *38:1* means reads resolve well to species level.
  #let w = 260pt
  #set align(center)
  #grid(columns: (1fr, 1fr), gutter: 10pt,
  // [
  //   #figure(
  //     image("../../out/figs/family_abundance_bar.png", width: w),
  //     caption: [\ Family log-abundances.]
  //   ) <family_ab>
  // ],
  [
    #figure(
      image("../../out/figs/genus_abundance_bar.png", width: w),
      caption: [Genus log-abundances.]
    ) <genus_ab>
  ],[
    #figure(
      image("../../out/figs/species_abundance_bar.png", width: w),
      caption: [Top 90% Species log-abundances.]
    ) <species_ab>
  ])
  
] 

#slide(title:"Results: Edit Distance MDS")[
  #set text(.9em)
  #grid(columns: (1fr, 1.5fr), gutter: 0pt)[
    - Pairwise Levenshtein distances on stratified 10,000-read sample
    - PERMANOVA
      - pseudo-F = 133.5, *p = .001*
        - Species membership is a strong predictor of read-sequence similarity
    - MDS reveals clear clustering by species
  ][
    #figure(
      align(right,image("../../out/figs/top_perc_read_str_dist_mds.png", width: 90%)),
      caption: [MDS of read Levenshtein distances, colored by species.]
    ) <mds>
  ]
]

#slide(title:"Discussion")[
  - SARS-CoV-2 and SARS CoV Tor2 rank among top 5 most abundant species (@species_ab)
    - replicates Zhou et al.: _R. malayanus_ hosts CoV @Zhou2020
    - Supports _Rhinolophidae_ are CoV reservoirs @Han2023
  - Virome dominated by bacteriophages
  - Edit distance analysis provides complementary confirmation
    - reads from the same species are more similar to each other than to reads from other species (@mds)
  - Shotgun metagenomics is a valid tool for monitoring pathogen zoonosis
]

#show bibliography: set text(.8em)
#let bib = bibliography("../bibliography.bib")
#bibliography-slide(bib, extra: align(bottom)[ https://github.com/canderson318/metagenome_profiler])
