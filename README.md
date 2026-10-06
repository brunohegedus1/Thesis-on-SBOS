# Thesis on SBOS

MSc thesis of Bruno Muzy Hegedüs, "Implementation of Laser-Speckled
Background Oriented Schlieren to Wind Tunnel Experiments". The repository
holds the LaTeX source of the report, its figures, and the MATLAB code that
produces the results. The sections on the template, further down, come with
the TU Delft class the report is built on.

## Report source

Chapters sit in `mainmatter/`, appendices in `appendix/`, the front matter in
`frontmatter/` and every figure in `figures/`. `report.tex` inputs them all
and `tudelft-report.cls` is the document class. Build with
`latexmk -xelatex report.tex`.

## MATLAB implementations

`Final Implementations/` holds the MATLAB code behind the results, one folder
per experiment:

| Folder | What it produces |
|---|---|
| `Methodology` | the design space exploration of Chapter 4 and the surface maps of Appendix C |
| `Angled glass experiment` | sensitivity against defocus distance, Figure 6.1 |
| `Speckle Pattern investigation` | speckle size per point and its uncertainty, Tables 6.7 and 6.8 |
| `Compressible jet experiment` | displacement fields, density and density gradient fields, shock cell lengths |
| `Aerospike experiment` | density gradient fields and speckle size of the aerospike runs |

`Final Implementations/README.md` lists every script, what it does, the data
it reads and what was left out of the collection.

### Measurement data

The DaVis exports are not in this repository. They run from 30 to 160 MB per
file and 7.4 GB in total, past the 100 MB per file that GitHub accepts, so
`.gitignore` keeps these out and they live on the local disk only:

* `Compressible jet experiment/CC results/` and `BOS_data/`
* `Aerospike experiment/Aerospike Data/` and `Aerospike Speckle Data/`
* `Speckle Pattern investigation/I*.csv`, `I*_2.txt`, `I3_2.mat` and `Exp1/`

What is tracked is every script plus the small files the plotting scripts
read: the sweep outputs of Chapter 4, the angled glass measurements, the
speckle size tables, the digitised Emden and Panda curves and the extracted
centrelines, 0.8 MB together. Those scripts therefore run from a fresh clone.
The ones that read a raw export need the data put back next to them, at the
paths listed above.

## Report template

This template aims
 to simplify and improve the (Xe)LaTeX report/thesis template by Delft University of Technology with the following three main design principles:

* **Simplicity First:** A class file that has been reduced by nearly 70% to simplify customization;
* **Effortless:** A careful selection of common packages to get started immediately;
* **Complete:** Ready-to-go when it comes to the document and file structure.

This template works with _pdfLaTeX_, _XeLaTeX_ and _LuaLaTeX_. In order to adhere to the TU Delft house style, either _XeLaTeX_ or _LuaLaTeX_ is required, as it supports TrueType and OpenType fonts. _BibLaTeX_ is used for the bibliography with as backend _biber_. Please visit https://dzwaneveld.github.io/report/ for the full documentation.

Cover | Title | Chapter
--- | --- | ---
<img src="https://dzwaneveld.github.io/images/report-template.jpg"> | <img src="https://dzwaneveld.github.io/images/report-template-title.jpg"> | <img src="https://dzwaneveld.github.io/images/report-template-chapter.jpg">

## Documentation (Abridged)

As a report/thesis is generally a substantial document, the chapters and appendices have been separated into different files and folders for convenience. The folders are based on the three parts in the document: the frontmatter, mainmatter and appendix. All files are inserted in the main file, `report.tex`, using the `\input{filename}` command. The document class, which can be found in `tudelft-report.cls`, is based on the `book` class.

The template will automatically generate a cover when the `\makecover` command is used. The title, subtitle and author will also be present on the title page. To give greater flexibility over the title page, the layout is specified in `title-report.tex`. A title page for theses is also available: `title-thesis.tex`. Change the corresponding `\input{...}` command in the main file to switch. 

The bibliography has been set up in `report.tex` to allow for easy customization. It is included in the table of contents and renamed to 'References' using the `heading=bibintoc` and `title=References` options of the `\printbibliography` command respectively. If you would like to use a different `.bib` file, change the command `\addbibresource{report.bib}` accordingly. 

*→ Visit https://dzwaneveld.github.io/report/ for the full documentation.*

## License

This [report/thesis template](https://github.com/dzwaneveld/TU-Delft-Unofficial-Report-Template) by Daan Zwaneveld is licensed under [CC BY-NC 4.0](https://creativecommons.org/licenses/by-nc/4.0/). No attribution is required in PDF outputs created using this template.
