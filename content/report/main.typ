
#import "utils.typ": *
#import "@preview/zebraw:0.6.1": *
#show: zebraw

////////////////////
// Document Settings
////////////////////
#set page(margin: .5in, width: 8.5in, height: 11in, numbering: "1")
#let fntsz = 11pt
#let after = fntsz * 1.5
#set text(font: "Georgia", size: fntsz)
#set par(leading: fntsz, spacing: after)
// #set heading(numbering: "1.")

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
    set text(style:"italic")
    it
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

//////////////////////
//////////////////////
//////// BODY ////////
//////////////////////
//////////////////////

// Written report (25%): This maximum 4-page (11pt arial font, ½ inch margins) technical document with an introduction, methods, and results sections. The introduction gives any necessary background and describes and justifies your strategy. The method section defines the input, presents the algorithm, and details the expected output. The results contain experiments, results, and discussion. You should cite the appropriate scientific literature where appropriate. The list of citations does not count against the page limit.

#title[Metatranscriptomic Virome Characterization of _Rhinolophus malayanus_ Bat Feces]
#align(center)[By: Christian Anderson\ #today.display("[day] [month repr:long] [year]")]


#align(bottom, outline())
#pagebreak()


= Introduction
// - assigned problem
// - Overview 
// - Computational translation
// - Biological Relevance
// - importance 
// - example 

Environmental monitoring is a critical component of pandemic readiness. The objective of the programming assignment is to develop 


- covid pandemic
- importance of environmental monitoring
- introduction of metagenome sampling
    - shotgun sequencing
    - host + microbiome
- profile microbiome by removing host signal
- look at community abundances
    - what communities are present
    - which dominate
    - (why do abundances matter?)
- need simple method for characterizing virome
- example with bats
    - cite @Zhou2020

= Methods
I developed a general method for characterizing the viral communities in a metagenomic sample where the input is a short-read fastq file containing short-read sequences and the output is a characterization of the virome communities present in the sample. 

- host alignment and read filtering
- phred filtering
- viral database querying 
- abundance analysis
- read edit distance analysis

= Results

- how many original reads
- how many removed because 
    - from host
    - low phred
- how many reads were classified
- rank classificationo
    - more species than any other so reads are distinct enough for accurate LCA 
- read edit distance
    - permanova
    - reads cluster with species
        - species is main driver of read patterns

= Discussion
- 

= AI use statement
AI was used in the following cases during code development.

- Helped me understand file IO byte offsets and how to find specific line's byte offsets using ```python open(file).tell()```
- Helped me understand how the line pointer shifts each time ```python open(file).readline()``` is called
- Helped me understand how to use an initializer with ```python Multiprocessing.Pool()``` to share an object globally across workers
- 



#pagebreak()
#bibliography("../bibliography.bib", title: "References")