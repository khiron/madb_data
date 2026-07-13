- Cover page
    [] change date to 2025
    
- Introduction
- Methods
- Results
- Discussion
- Bibiography
  - more citations
  - fix de Bruijn 
  - fix WELCH
- Figures & Tables
    - General
        - First sentence in caption should be a statement that tells the reader what they're meant to take away from the data
        - Don't use bold in figure captions
        - Don't use red and green
```python
species_colors = {
    "Chimpanzee":    "#CC0000",   # Bold Red (still legible for most CVD types)
    "Gorilla":    "#E69F00",  # Orange
    "Macaque": "#0072B2",  # Blue
}
```
```latex
\usepackage{xcolor}

% Define species colors
\definecolor{chimpanzeeRed}{HTML}{CC0000}
\definecolor{gorillaOrange}{HTML}{E69F00}
\definecolor{macaqueBlue}{HTML}{0072B2}
```
        - figures same size/font/colours

    - Figure 1.1 [debruijngraph_example.png]

    - Figure 1.2 [needleman_wunsch_1.png,needleman_wunsch_2.png]
Don't make the caption bold-faced. You need to make the individual figures or subplots consistent between these two figures. You need to describe what the arrows mean and what the individual values in the cells mean

    - Figure 1.3 [smith_waterman.png]

Don't make the caption bold-faced. You need to make the individual figures or subplots consistent between these two figures. You need to describe what the arrows mean and what the individual values in the cells mean    

    - Figure 2.1 [pdist_pairs.png]

        - Turn into a pdf
        - use solution for MathJax 
        - not bold-faced
        - add n of each pair
        - identify variabls using $var$
        - update colours
        - mention that the n's are teh same because this is from a sample that contain all sequences

    - Figure 2.2
    - Figure 2.3
    - Figure 2.4
    - Figure 2.5
    - Figure 2.6
    - Table 2.2 
    - Table 2.3
        - add n as a column
