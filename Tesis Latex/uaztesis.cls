% uaztesis.cls
% Archivo de estilo LaTeX para formato de tesis de la 
% Unidad Académica de Ingeniería Eléctrica de la UAZ, 
% adaptado por Manuel Reta Hernández de los formatos 
% de tesis de UW-Madison (Eric Benedict) y de the Purdue 
% University (Dave Kraynie / Darrell McCauley)
%
% Abril, 2005
%=============================================================================
%\NeedsTeXFormat{LaTeX2e}
%\ProvidesClass{uaztesis}[abril de 2005, Univ. Aut. de Zacatecas]
%\typeout{Document Class `uaztesis' }

% Uso de paquetes 

\usepackage[latin1]{inputenc}
\usepackage{amssymb,amsmath,amsfonts}
\usepackage{color,graphicx,curves}
\usepackage{epsfig}
\usepackage{times}

\renewcommand{\ttdefault}{cmtt}    

%=============================================================================
% Tamaño de papel
%=============================================================================
\setlength\paperheight{11in}
\setlength\paperwidth{8.5in}

%=============================================================================
% La opción de tamaño de letra es manejada leyendo 
% un archivo diferente para cada tamaño, como sigue,
% las instrucciones para definir el tamaño de letra:
% 10pt : withe10,    11pt : withe11,    12pt : withe12
% Ejecutadas definiendo \@ptsize para manejar el último dígito
%=============================================================================
\newcommand\@ptsize{}
\DeclareOption{10pt}{\renewcommand\@ptsize{0}}  % opción 10pt lee en withe10.STY
\DeclareOption{11pt}{\renewcommand\@ptsize{2}   
     \typeout{11pt not available, subsituting 12pt}}  % opción 11pt lee en withe12.STY
\DeclareOption{12pt}{\renewcommand\@ptsize{2}}        % opción 12pt lee en withe12.STY

%=============================================================================
% Impresión a uno o dos lados de la página 
% \@twosidefalse                      % Default es a una página.
%=============================================================================
\DeclareOption{twoside}{\@twosidetrue % Define opción a ambos lados
           \@mparswitchtrue}          % Marginpars rebazan los márgenes de la página

%=============================================================================
%  Opción de checar márgenes 
% \overfullrule = 0pt                 % Default es no marcar extralímite de hboxes
%=============================================================================
\DeclareOption{margincheck}{\overfullrule 5pt}  % Marca extralímite de hboxes

%=============================================================================
%  Opción de margen de tesis
%  \@msmarginfalse                    % Default es márgenes de 1"
%  Los márgenes deben ser modificados en \documentstyle ya que se 
%  especifican they cuando se procesa el prámbulo.  La bandera \@msmargin
%  también se usa para cabiar los defaults de \thesis y \degree
%=============================================================================
\newif\if@msmargin                    % Crea un switch, default es false
\DeclareOption{tesis}{\@msmargintrue  % Define márgen izquierdo de 1.5"
      \typeout{Usando márgenes y defaults de Tesis de maestría}}

%=============================================================================
% La instrucción \ProcessOptions\relax hace que se ejecute cada 
% instrucción definida con \DeclareOption{FOO} y para la cual el usuario 
% tecleo la opción FOO en su \documentclass. Por cada opción BAR tecleada 
% para lo cual BAR no está definida, el archivo BAR.sty será leido después 
% de que se ejecute el presente archivo uaztesis.sty 
%=============================================================================
\ExecuteOptions{10pt}
\ProcessOptions\relax
\RequirePackage{withe1\@ptsize}

%=============================================================================
%  Cambia la definición de fuentes a la estandard LaTeX
%=============================================================================
\DeclareOldFontCommand{\rm}{\normalfont\rmfamily}{\mathrm}
\DeclareOldFontCommand{\sf}{\normalfont\sffamily}{\mathsf}
\DeclareOldFontCommand{\tt}{\normalfont\ttfamily}{\mathtt}
\DeclareOldFontCommand{\bf}{\normalfont\bfseries}{\mathbf}
\DeclareOldFontCommand{\it}{\normalfont\itshape}{\mathit}
\DeclareOldFontCommand{\sl}{\normalfont\slshape}{\@nomath\sl}
\DeclareOldFontCommand{\sc}{\normalfont\scshape}{\@nomath\sc}

%=============================================================================
% LISTAS
%=============================================================================
%  ENUMERATE
%  La numeración se hace con 4 contadores: enumi, enumii, enumiii
%  y enumiv, donde enumN controla la numeración del 
%  nivel Nth. La etiqueta es generada por las instrucciones 
%  \labelenumi... \labelenumiv.  La expansión de 
%  \p@enumN\theenumN define la salida de una instrucción \ref
%=============================================================================
\def\labelenumi{\arabic{enumi}.}
\def\theenumi{\arabic{enumi}}
\def\labelenumii{(\alph{enumii})}
\def\theenumii{\alph{enumii}}
\def\p@enumii{\theenumi}
\def\labelenumiii{\roman{enumiii}.}
\def\theenumiii{\roman{enumiii}}
\def\p@enumiii{\theenumi(\theenumii)}
\def\labelenumiv{\Alph{enumiv}.}
\def\theenumiv{\Alph{enumiv}}
\def\p@enumiv{\p@enumiii\theenumiii}

%=============================================================================
% ITEMIZE
% La lista con símbolos es controlada por 4 instrucciones:
% \labelitemi, \labelitemii, \labelitemiii, y \labelitemiv, 
% Las cuales definen las etiquetas de varios niveles. 
%=============================================================================
\def\labelitemi{$\bullet$}
\def\labelitemii{\bface --}
\def\labelitemiii{$\ast$}
\def\labelitemiv{$\cdot$}

%=============================================================================
%  VERSE
%  El ambiente verse se define haciendo uso de los parámetros
%  de ambiente. Se teclea \\ para finalizar una línea.
%  Esto se implementa con \let en \\ equal \@centercr.
%=============================================================================
\def\verse{\let\\=\@centercr
  \list{}{\itemsep\z@ \itemindent -1.5em\listparindent \itemindent
          \rightmargin\leftmargin\advance\leftmargin 1.5em}\item[]}
\let\endverse=\endlist

%=============================================================================
%  QUOTATION
%  Llena lineas; indenta párrafos
%=============================================================================
\def\quotation{\list{}{\listparindent 1.5em
  \itemindent\listparindent
  \rightmargin\leftmargin\parsep 0pt plus 1pt}\item[]}
\let\endquotation=\endlist

%=============================================================================
%  QUOTE - igual que quotation excepto que no indenta párrafos
%=============================================================================
\def\quote{\list{}{\rightmargin\leftmargin}\item[]}
\let\endquote=\endlist

%=============================================================================
%  DESCRIPTION
%  Para cambiar el formato de la etiqueta, redefinir \descriptionlabel.
%=============================================================================
\def\descriptionlabel#1{\hspace\labelsep \bface #1}
\def\description{\list{}{\labelwidth\z@ \itemindent-\leftmargin
  \let\makelabel\descriptionlabel}}
\let\enddescription\endlist
\newdimen\descriptionmargin
\descriptionmargin=3em

%=============================================================================
% OTROS AMBIENTES
%=============================================================================
% THEOREM
% \@begintheorem ... \@endtheorem son las instrucciones 
% ejecutadas al inicio y fin de un ambiente tipo 
% teorema (definido por el usuario).
% Except \@opargbegintheorem is executed when an optional argument is
% given.  Cf. LATEX.TEX.
%=============================================================================
\def\@begintheorem#1#2{\trivlist \item[\hskip \labelsep{\bface #1\ #2}]}
\def\@opargbegintheorem#1#2#3{\trivlist
  \item[\hskip \labelsep{\bface #1\ #2\ (#3)}]}
\def\@endtheorem{\endtrivlist}

%=============================================================================
%  PÁGINA DE PRESENTACIóN
%  En los ambientes normales, el ambiente del título de 
%  página sólo comienza y termina una página, inhibiendo 
%  los números de página. También reinicia el número 
%  de página a cero. Aún en estilo de dos columnas hace una 
%  página de una columna. 
%=============================================================================
\newenvironment{titlepage}
 {\@restonecolfalse\if@twocolumn\@restonecoltrue\onecolumn
  \else \newpage \fi \thispagestyle{empty}
% \c@page\z@ -- borrada la cuenta de páginas
}{\if@restonecol\twocolumn \else \newpage \fi}

%=============================================================================
% ARRAY Y TABULAR
%=============================================================================
\arraycolsep 5pt     % Half the space between columns in an array environment.
\tabcolsep 6pt       % Half the space between columns in a tabular environment.
\arrayrulewidth .4pt % Width of rules in array and tabular environment.
\doublerulesep 2pt   % Space between adjacent rules in array or tabular env.

%=============================================================================
% TABBING
%=============================================================================
\tabbingsep \labelsep   % Espacio usado por la 
                        % instrucción \' (Ver manual LaTeX)

%=============================================================================
%  MINIPAGE
%  \@minipagerestore is called upon entry to a minipage environment to
%  set up things that are to be handled differently inside a minipage
%  environment. In the current styles, it does nothing.
%
%  \skip\@mpfootins : plays same role for footnotes in a minipage as
%                    \skip\footins does for ordinary footnotes
%=============================================================================
\skip\@mpfootins = \skip\footins

%=============================================================================
% FRAMEBOX
%=============================================================================
\fboxsep = 3pt    % Espacio entre caja y texto por \fbox y \framebox.
\fboxrule = .4pt  % Ancho de reglas en caja hecho por \fbox y \framebox.

%=============================================================================
% CAPÍTULOS Y SECCIONES
%=============================================================================
% Definiendo contadores:
% \newcounter{NEWCTR}[OLDCTR] : Define NEWCTR como un contador, el cual 
%                               es reiniciado en cero cuando el contador 
%                               OLDCTR es brincado.
%                               El contador OLDCTR debió haberse definido
%=============================================================================
\newcounter {part}
\newcounter {chapter}
\newcounter {section}[chapter]
\newcounter {subsection}[section]
\newcounter {subsubsection}[subsection]
\newcounter {paragraph}[subsubsection]
\newcounter {subparagraph}[paragraph]
\newcounter {tocpage}
\newcounter {lofpage}
\newcounter {lotpage}
\newcounter {listofheading}

%=============================================================================
% Para cualquier contador CTR, \theCTR es un macro que define la versión impresa
% del contador CTR.  Está definido en términos de los siguientes macros:
%
%  \arabic{COUNTER} : El valor de COUNTER se imprime como número arábigo.
%  \roman{COUNTER}  : Se imprime como un número romano en minúsculas.
%  \Roman{COUNTER}  : Se imprime como un número romano en mayúsculas.
%  \alph{COUNTER}   : Se imprime como una letra minúscula:
%                         1 = a, 2 = b, etc.
%  \Alph{COUNTER}   : Se imprime como una letra mayúscula
%                         1 = A, 2 = B, etc.
%=============================================================================
\def\thepart          {\Roman{part}}
\def\thechapter       {\arabic{chapter}}
\def\thesection       {\thechapter.\arabic{section}}
\def\thesubsection    {\thesection.\arabic{subsection}}
\def\thesubsubsection {\thesubsection.\arabic{subsubsection}}
\def\theparagraph     {\thesubsubsection.\arabic{paragraph}}
\def\thesubparagraph  {\theparagraph.\arabic{subparagraph}}

%=============================================================================
% TABLA DE CONTENIDOS, ETC.
%=============================================================================
% Una instrucción \subsection escribe una instrucción
%       \contentsline{subsection}{TITLE}{PAGE} 
% en el archivo .toc, donde TITLE indica el contenido 
% y PAGE indica el número de página. Si las subsecciones
% son numeradas, entonces TITLE será de la forma 
%       \numberline{NUM}{HEADING}
% donde NUM es el número producido por \thesubsection.  
% Otras instrucciones de seccionado trabajan en forma similar.
%
% Una instrucción \caption en un ambiente de 'figure' escribe
%   \contentsline{figure}{\numberline{NUM}{CAPTION}}{PAGE}
% en el archivo .lof, donde NUM es el número producido por
% \thefigure y CAPTION es el nombre de la figura.  
% Trabaja en forma similar para el ambiente 'table'
%
% La instrucción \contentsline{NAME} expande a \l@NAME.  
% Por tanto, para especificar el contenido general se debe 
% definir \l@chapter, \l@section, \l@subsection, ... 
% Para especificar la lista de figuras, se debe definir 
% \l@figure y así sucesivamente. La mayoría de ellos se puede 
% definir con la instrucción \@dottedtocline, la cual 
% trabaja como sigue:
%
% \@dottedtocline{LEVEL}{INDENT}{NUMWIDTH}{TITLE}{PAGE}
%    LEVEL    : An entry is produced only if LEVEL < or = value of
%               'tocdepth' counter.  Note, \chapter is level 0, \section
%               is level 1, etc.
%    INDENT   : The indentation from the outer left margin of the start of
%               the contents line.
%    NUMWIDTH : The width of a box in which the section number is to go,
%               if TITLE includes a \numberline command.
%
% Esta instrucción usa los siguientes parámetros, los cuales son puestos
% con \def (por tanto se pueden usar em's para hacerlos depender del tipo de letra).
%   \@pnumwidth : ancho de una caja en la cual se pone el número de página.
%   \@tocrmarg  : márgen derecho para múltiples entradas. 
%                 Uno quiere \@tocrmarg > ó = \@pnumwidth
%   \@dotsep    : Separación entre puntos, en unidades mu. Debe ser 
%                 \def'd para un número como 2 ó 1.7
%=============================================================================
\def\@pnumwidth{1.55em}
\def\@tocrmarg {2.55em}
\def\@dotsep{4.5}
\setcounter{tocdepth}{2}

%=============================================================================
% The value of \@tocchapterskip when the table of contents is
% formatted is the only value that matters.  Changing this value on
% the fly probably won't do what you want.  Re-defining \l@chapter has
% the same problem.  The only way around the problem I've thought of
% is to write the spacing commands directly into the toc file from the
% chapter command, as is done by \appendix, for instance.
%=============================================================================
\def\@tocchapterskip{1.0em plus 1pt}

\def\l@part#1#2{\addpenalty{-\@highpenalty}
  \addvspace{2.25em plus 1pt}  % space above part line
  \begingroup
  \@tempdima 3em         % width of box holding part number, used by
    \parindent \z@ \rightskip \@pnumwidth             %% \numberline
    \parfillskip -\@pnumwidth
    {\large \bface       % set line in \large boldface
    \leavevmode          % TeX command to enter horizontal mode.
    #1\hfil \hbox to\@pnumwidth{\hss #2}}\par
    \nobreak             % Never break after part entry
  \endgroup}

\def\l@chapter#1#2{\pagebreak[3]
  \vskip\@tocchapterskip % space above chapter line
  \@dottedtocline{0}{0pt}{1.5em}{{\bface #1}}{#2}}

\def\l@section{\@dottedtocline{1}{1.5em}{2.3em}}
\def\l@subsection{\@dottedtocline{2}{3.8em}{3.2em}}
\def\l@subsubsection{\@dottedtocline{3}{7.0em}{4.1em}}
\def\l@paragraph{\@dottedtocline{4}{10em}{5em}}
\def\l@subparagraph{\@dottedtocline{5}{12em}{6em}}
\def\l@figure{\vskip 1.0em plus 1pt \@dottedtocline{0}{0em}{2.5em}}
\let\l@table\l@figure

%=============================================================================
% TABLA DE CONTENIDOS; LISTA DE FIGURAS; LISTA DE TABLAS
%=============================================================================
% In report style, \tableofcontents, \listoffigures, etc. are always
% set in single-column style.  @restonecol is used to keep track of
% whether we need to switch back to double column style after the toc.
%
% The only known problem now is that the first page with the new
% layout is too long.  The problem seems to be that the change to
% textheight doesn't take place on the first page.  Even if it's the
% first line in the table of contents macro.  Presumably the same
% problem also occurs in the lof and lot.
%
% I'm taking a shot at fixing the problem by dropping in a throw-away
% page between the change to the height parameters and the start of
% the chapter.  Isn't elegance wonderful?
%
%=============================================================================
\def\tableofcontents{\@tableof{Contenido General}{}{tocpage}{toc}{}}

\def\listoffigures{
  \@tableof{Lista de figuras}{Figura}{lofpage}{lof}
  {\protect\addcontentsline{toc}{chapter}{Lista de figuras}}}

\def\listoftables{
  \@tableof{Lista de tablas}{Tabla}{lotpage}{lot}
  {\protect\addcontentsline{toc}{chapter}{Lista de tablas}}}

\def\@tableof#1#2#3#4#5{
{ % limit scope of following declarations!!
  \@restonecolfalse\if@twocolumn\@restonecoltrue\onecolumn\fi
  \addtolength{\textheight}{-30pt}       % -24-16
  \addtolength{\majorheadskip}{-30pt}    % -24-16
  \addtolength{\headheight}{25 pt}       %  36+16   
  \addtolength{\headsep}{-15pt}          % -12   
%  \separatorpage{DISCARD THIS PAGE}
  \chapter*{#1}
  #5
  \relax\markboth{#1}{#1}
  \hbox to \hsize{#2 \hfil Pag.}
  \singlespace
  \setcounter{#3}{0}
  \setcounter{listofheading}{1}  % change from 0 to 1 by mccauley, 14may93
  \def\@oddhead{\vbox to \headheight{\vspace{4pt}
    \hbox to \hsize{\hfil\rm\thepage} \vfil
    \ifnum\value{#3}=1
%      \ifnum\value{listofheading}=2
%        \hbox to \hsize{Apéndice\hfil} \vspace{4pt} \fi
      \ifnum\value{listofheading}=1
        \stepcounter{listofheading} \fi
      \hbox to \hsize{#2 \hfil Pag.}
    \else
      \setcounter{#3}{1}
    \fi}}
  \def\@evenhead{\vbox to \headheight{\vspace{4pt}
    \hbox to \hsize{\rm\thepage\hfil} \vfil
    \ifnum\value{#3}=1
%      \ifnum\value{listofheading}=2
%        \hbox to \hsize{Apéndice\hfil} \vspace{4pt} \fi
      \ifnum\value{listofheading}=1
        \stepcounter{listofheading} \fi
      \hbox to \hsize{#2 \hfil Pag.}
    \else
      \setcounter{#3}{1}
    \fi}}
  \@starttoc{#4}  \if@restonecol\twocolumn\fi
  \newpage
}}

%=============================================================================
% REFERENCIAS
%=============================================================================
% The thebibliography environment executes the following commands:
%
%  o start a new 'chapter' with BIBLIOGRAPHY as the heading
%  o produce a separator page for the bibliography
%
%  \def\newblock{\hskip .11em plus .33em minus -.07em} --
%      Defines the `closed' format, where the blocks (major units of
%      information) of an entry run together.
%
%  \sloppy  -- Used because it's rather hard to do line breaks in
%      bibliographies,
%
%  \sfcode`\.=1000\relax --
%      Causes a `.' (period) not to produce an end-of-sentence space.
%=============================================================================
% \altbibtitle
%   The default title for the References chapter is ``LIST OF REFERENCES''
%   Since some people prefer ``BIBLIOGRAPHY'', the command
%   \altbibtitle has been added to change the chapter title.
%   This command does nothing more than change REFERENCES to BIBLIOGRAPHY
%============================================================================
\def\@bibchaptitle{Referencias}
%\def\altbibtitle{\def\@bibchaptitle{Bibliografía}}
\def\thebibliography#1{
  %\separatorpage{\@bibchaptitle}
  \global\@bibpresenttrue
  \chapter*{\@bibchaptitle\markboth{\@bibchaptitle}{\@bibchaptitle}}
  \addcontentsline{toc}{chapter}{\@bibchaptitle}
  \vspace{0.375in}    % added to match 4 line requirement
  \interlinepenalty=10000 % added to prevent breaking of bib entries
  \singlespace\list
  {[\arabic{enumi}]}{\settowidth\labelwidth{[#1]}\leftmargin\labelwidth
    \advance\leftmargin\labelsep \usecounter{enumi}}
  \def\newblock{\hskip .11em plus .33em minus -.07em}
  \sloppy
  \sfcode`\.=1000\relax}
\let\endthebibliography=\endlist

%=============================================================================
% \def\@biblabel#1{[#1]\hfill}  % Produces the label for \bibitem[...]
% \def\@cite#1{[#1]}            % Produces the output of the \cite command.
%=============================================================================

%=============================================================================
%   EL AMBIENTE ÍNDICE
%   Produces double column format, with each paragraph a separate entry.
%   The user commands \item, \subitem and \subsubitem are used to
%   produce the entries, and \indexspace adds an extra vertical space
%   that's the right size to put above the first entry with a new letter
%   of the alphabet.
%=============================================================================
\newif\if@restonecol
\def\theindex{\@restonecoltrue\if@twocolumn\@restonecolfalse\fi
  \columnseprule \z@
  \columnsep 35pt\twocolumn[\@makeschapterhead{Index}]
  \markboth{Índice}{Índice}\thispagestyle{plain}\parindent\z@
  \parskip\z@ plus .3pt\relax\let\item\@idxitem}
\def\@idxitem{\par\hangindent 40pt}
\def\subitem{\par\hangindent 40pt \hspace*{20pt}}
\def\subsubitem{\par\hangindent 40pt \hspace*{30pt}}
\def\endtheindex{\if@restonecol\onecolumn\else\clearpage\fi}
\def\indexspace{\par \vskip 10pt plus 5pt minus 3pt\relax}

%=============================================================================
% NOTAS DE PIE DE PÁGINA
%=============================================================================
% \footnoterule is a macro to draw the rule separating the footnotes from
% the text.  It should take zero vertical space, so it needs a negative
% skip to compensate for any positive space taken by the rule.  (See
% PLAIN.TEX.)
%=============================================================================
\def\footnoterule{\kern-3\p@
  \hrule width .4\columnwidth
  \kern 2.6\p@}                 % The \hrule has default height of .4pt .
% \newcounter{footnote}
\@addtoreset{footnote}{chapter}  % Numbers footnotes within chapters

%=============================================================================
% \@makefntext{NOTE} :
%   Must produce the actual footnote, using \@thefnmark as the mark
%   of the footnote and NOTE as the text.  It is called when effectively
%   inside a \parbox of width \columnwidth (i.e., with \hsize =
%   \columnwidth).
%
%   The following macro indents all lines of the footnote by 10pt, and
%   indents the first line of a new paragraph by 1em.  To change these
%   dimensions, just substitute the desired value for '10pt' [in both
%   places] or '1em'.  The mark is flushright against the footnote.
%     \long\def\@makefntext#1{\@setpar{\@@par\@tempdima \hsize
%        \advance\@tempdima-10pt\parshape \@ne 10pt \@tempdima}\par
%        \parindent 1em\noindent \hbox to \z@{\hss$^{\@thefnmark}$}#1}
%
%   A simpler macro is used, in which the footnote text is
%   set like an ordinary text paragraph, with no indentation except
%   on the first line of a paragraph, and the first line of the
%   footnote.  Thus, all the macro must do is set \parindent
%   to the appropriate value for succeeding paragraphs and put the
%   proper indentation before mark.
%=============================================================================
\long\def\@makefntext#1{\parindent 1em\noindent
  \hbox to 1.8em{\hss$^{\@thefnmark}$}\singlespace\footnotesize#1}

%=============================================================================
% \@makefnmark :
%   Macro to generate the footnote marker that goes in the text. Default used.
%=============================================================================

%=============================================================================
% FIGURAS Y TABLAS
%=============================================================================
% Float placement parameters.  See LaTeX manual for their definition.
%
% Note: some really flaky hiccups can happen because of interactions between
% these numbers.  For instance, if \textfraction and \floatpagefraction
% don't have any overlap figures can be drawn that are too large to fit
% on a text page and too small to fit on a float page.  This problem causes
% all of the figures to go to the end of the chapter.
%=============================================================================
\setcounter{topnumber}{2}
\def\topfraction{.5}
\setcounter{bottomnumber}{1}
\def\bottomfraction{.3}
\setcounter{totalnumber}{3}
\def\textfraction{.5}     % previously .55
\def\floatpagefraction{.3}
\setcounter{dbltopnumber}{2}
\def\dbltopfraction{.7}
\def\dblfloatpagefraction{.5}

%=============================================================================
% \@makecaption{NUMBER}{TEXT} : Macro to make a figure or table caption.
%      NUMBER : Figure or table number--e.g., 'Figure 3.2'
%      TEXT   : The caption text.
%  Macro should be called inside a \parbox of right width, with \normalsize.
%=============================================================================
\long\def\@makecaption#1#2{
  \vskip 10pt
  \singlespace\centering
  #1\ \ #2}

%=============================================================================
% ORIGINAL \@makecaption - Centering OR Flushleft if long.
% Thesis office doesn't really like this version.
%=============================================================================
% \long\def\@makecaption#1#2{
%   \vskip 10pt
%   \setbox\@tempboxa\hbox{#1\ \ #2}
%   \ifdim \wd\@tempboxa >\hsize           % IF   longer than one line:
%     {\singlespace \unhbox\@tempboxa\par} % THEN set as ordinary paragraph.
%   \else                                  % ELSE center.
%        \hbox to\hsize{\hfil\box\@tempboxa\hfil}
%   \fi}
%=============================================================================

%=============================================================================
% NEW \@makecaption - Thesis Format may not like the new format
%=============================================================================
% \newlength\@capnumlen
% \newlength\@captionlen
% \long\def\@makecaption#1#2{
%   \vskip 10pt
%   \setbox\@tempboxa\hbox{#1\ \ #2}
%   \ifdim \wd\@tempboxa >\hsize           % IF   longer than one line:
%   { \parindent 0pt
%     \singlespace\normalsize
%     \setlength{\@captionlen}{\hsize}
%     \settowidth{\@capnumlen}{#1\ \ }
%     \addtolength{\@captionlen}{-\@capnumlen}
%     #1\ \ \parbox[t]{\@captionlen}{\noindent#2}
%     \par\nobreak}                        % TeX penalty to prevent page break.
%   \else                                  % ELSE center.
%        \hbox to\hsize{\hfil\box\@tempboxa\hfil}
%   \fi}
%=============================================================================

%=============================================================================
% To define a float of type TYPE (e.g., TYPE = figure), the document style
% must define the following.
%
%  \fps@TYPE   : The default placement specifier for floats of type TYPE.
%
%  \ftype@TYPE : The type number for floats of type TYPE.  Each TYPE has
%                associated a unique positive TYPE NUMBER, which is a power
%                of two.  E.g., figures might have type number 1, tables type
%                number 2, programs type number 4, etc.
%
%  \ext@TYPE   : The file extension indicating the file on which the
%                contents list for float type TYPE is stored.  For example,
%                \ext@figure = 'lof'.
%
%  \fnum@TYPE  : A macro to generate the figure number for a caption.
%                For example, \fnum@TYPE == Figure \thefigure.
%
%  The actual float-making environment commands--e.g., the commands
%  \figure and \endfigure--are defined in terms of the macros \@float
%  and \end@float, which are described below.
%
%  \@float{TYPE}[PLACEMENT] : Macro to begin a float environment for a
%     single-column float of type TYPE with PLACEMENT as the placement
%     specifier.  The default value of PLACEMENT is defined by \fps@TYPE.
%     The environment is ended by \end@float.
%     E.g., \figure == \@float{figure}, \endfigure == \end@float.
%=============================================================================

%=============================================================================
% FIGURE
%=============================================================================
\newcounter{figure}[chapter]
\def\thefigure{\thechapter.\@arabic\c@figure}
\def\fps@figure{tbp}
\def\ftype@figure{1}
\def\ext@figure{lof}
\def\fnum@figure{\small{Figura} \thefigure}
\def\figure{\@float{figure}}
\let\endfigure\end@float
\@namedef{figure*}{\@dblfloat{figure}}
\@namedef{endfigure*}{\end@dblfloat}

%=============================================================================
% TABLE
%=============================================================================
\newcounter{table}[chapter]
\def\thetable{\thechapter.\@arabic\c@table}
\def\fps@table{tbp}
\def\ftype@table{2}
\def\ext@table{lot}
\def\fnum@table{\small{Tabla} \thetable}
\def\table{\@float{table}}
\let\endtable\end@float
\@namedef{table*}{\@dblfloat{table}}
\@namedef{endtable*}{\end@dblfloat}

%=============================================================================
% PÁGINA DE PRESENTACIÓN
%
%   - define \title{} \author{} \date{}
%   - The default degree is ``Doctor of Philosophy''
%     Degree can be changed using the command \degree{}
%   - for a masters project report, specify \project
%   - for a preliminary report, specify \prelim
%   - The default department is ``Electrical Engineering''
%     The department can be changed using the command \department{}
%   - once the above are defined, use \maketitle to generate the titlepage
%=============================================================================
\def\advisortitle#1{\gdef\@advisortitle{#1}}
\def\advisorname#1{\gdef\@advisorname{#1}}
\gdef\@advisortitle{empty}
\gdef\@advisorname{empty}
\def\gradoasesor2#1{\gdef\@gradoasesor2{#1}}
\def\nombreasesor2#1{\gdef\@nombreasesor2{#1}}
\gdef\@gradoasesor2{empty}
\gdef\@nombreasesor2{empty}
\newcommand\@thesistitlemedskip{0.25in}
\newcommand\@thesistitlebigskip{0.6in}
\newcommand{\degree}[1]{\gdef\@degree{#1}}
\newcommand{\report}{\gdef\@doctype{Reporte de Trabajo}}
\newcommand{\thesis}{\gdef\@doctype{Tesis de Maestría}}
\newcommand{\thesism}{\gdef\@doctype{Tesis de Maestría}}
\newcommand{\department}[1]{\gdef\@department{(#1)}}

\if@msmargin    % Set Defaults based on Margin definition
   \gdef\@degree{Ingeniero Electricista}   % Default is MS
   \gdef\@doctype{Tesis de Maestría}   % Default is Thesis
\else
   \gdef\@degree{Maestro en Ingeniería \\
   con Especialidad en Sistemas de Potencia} %Default is PhD
   \gdef\@doctype{Tesis de Maestría}         %Default is dissertation
\fi

\gdef\@department{Unidad Académica de Ingeniería 
Eléctrica} % Default is Electrical Engineering

\newcommand{\maketitle}{%
  \begin{titlepage}
%-----------------------------------------------------------------------------
% -- The thesis office doesn't like thanks on title page.  Put it in
% -- the acknowledgments.  This is here so you don't have to change
% -- your titlepage when converting from report style. -> from Purdue, but I
%        left it here since it seems compatible with UW-Madison, Eric
%-----------------------------------------------------------------------------
    \def\thanks##1{\typeout{Warning: `thanks' deleted from thesis titlepage.}}
    \let\footnotesize\small \let\footnoterule\relax \setcounter{page}{1}
%    \vspace*{0.2in}
    \begin{center}
    {\Large\bf\expandafter\uppercase\expandafter{UNIVERSIDAD 
    AUTÓNOMA DE ZACATECAS}}\\
     \vspace*{12pt}
      \begin{figure}[h]
      \centering
      \includegraphics[width=1.7in]{logoUAZnuevo.eps}\\
      \end{figure}
     \vspace*{36pt}
      {\bf\expandafter\uppercase\expandafter{\@title}} \\  
     \vspace*{12pt}
      {\bf\@author} \\ 
     \vspace*{48pt}
      \@doctype\ \\
      presentada a la \@department \\
      de acuerdo a los requerimientos de la Universidad 
      para obtener el título de \\
     \vspace*{48pt}
      {\bf\expandafter\uppercase\expandafter{\@degree}} \\
%       \\[\@thesistitlemedskip]
     \vspace*{20pt}
      \small{Directores de tesis: \\ 
      \@advisortitle\ \@advisorname\ y 
      \@gradoasesor2\ \@nombreasesor2}\\
     \vspace*{20pt}     
      \small{UNIDAD ACADÉMICA DE INGENIERÍA ELÉCTRICA \\
      \vspace*{8pt}
      Zacatecas, Zac.,} \@date \\
    \end{center}
  \end{titlepage}
  \setcounter{footnote}{0}
  \setcounter{page}{1} %title page is NOT counted
  \let\thanks\relax
  \let\maketitle\relax \let\degree\relax \let\project\relax \let\prelim\relax
  \let\department\relax
%  \gdef\@thanks{}\gdef\@degree{}\gdef\@doctype{}
%  \gdef\@department{}
%  \gdef\@author{}\gdef\@title{}
}

%=============================================================================
% HOJA DE REVISION DE TEMA DE TESIS
%=============================================================================

\def\dirtitle#1{\gdef\@dirtitle{#1}}
\def\dirname#1{\gdef\@dirname{#1}}
%\def\primrevtitle#1{\gdef\@primrevtitle{#1}}
%\def\primrevname#1{\gdef\@primrevname{#1}}
%\def\segrevtitle#1{\gdef\@segrevtitle{#1}}
%\def\segrevname#1{\gdef\@segrevname{#1}}
\def\fechasol#1{\gdef\@fechasol{#1}}
\def\fecharev#1{\gdef\@fecharev{#1}}
\def\carrera#1{\gdef\@carrera{#1}}
\gdef\@primrevtitle{empty}
\gdef\@primrevname{empty}
\gdef\@segrevtitle{empty}
\gdef\@segrevname{empty}
\gdef\@fechasol{empty}
\gdef\@fecharev{empty}
\gdef\@carrera{empty}

\def\revision{
             \thispagestyle{empty}
%                  \addtocounter{page}{-1}
\begin{center}
{\bf\expandafter\uppercase\expandafter{Aprobación de 
tema de Tesis de Maestría}}\\
\vspace*{12pt}
      \begin{figure}[h!]
      \centering
      \includegraphics[width=1.4in]{escudo.eps}\\
      \end{figure}
\end{center}
%\vspace*{8pt}

\@author

PRESENTE
\vspace*{18pt}

De acuerdo a su solicitud de tema de \@doctype\ del Programa de
\@carrera, con fecha \@fechasol, se acuerda asignarle el tema titulado:
\vspace*{6pt}
\begin{center}
{\bf\expandafter\uppercase\expandafter{\@title}}
\end{center}
\vspace*{6pt}

Se nombran revisores de Tesis a los profesores \@advisortitle\ \@advisorname\ 
y \@gradoasesor2\ \@nombreasesor2, notificándole a usted que dispone de un plazo máximo de seis meses,
a partir de la presente fecha, para la conclusión del documento final debidamente revisado. 

\vspace*{48pt}

\begin{center}
Atentamente \\
Zacatecas, Zac., \@fecharev \\
%\vspace*{12pt}
\vspace*{28pt}                       
\underline{\makebox[2.8in][r]{}}\\                  
\@dirtitle \hspace{3pt} \@dirname \\
Director de la \@department
\end{center}
}

\def\endrevision{\vfill \hfill\par\newpage}

%=============================================================================
% HOJA DE AUTORIZACION DE IMPRESION
%=============================================================================

\def\dirtitle#1{\gdef\@dirtitle{#1}}
\def\dirname#1{\gdef\@dirname{#1}}
\def\fechaaut#1{\gdef\@fechaaut{#1}}
\gdef\@dirtitle{empty}
\gdef\@dirname{empty}
\gdef\@fechaaut{empty}
\def\autorizacion{
             \thispagestyle{empty}
%                  \addtocounter{page}{-1}
\begin{center}
{\bf\expandafter\uppercase\expandafter{Autorización de 
impresión de Tesis de Maestría}}\\
\vspace*{12pt}
      \begin{figure}[h]
      \centering
      \includegraphics[width=1.4in]{escudo.eps}\\
      \end{figure}
\end{center}
%\vspace*{8pt}

\@author

PRESENTE
\vspace*{18pt}

%\vspace{24pt}
La Dirección de la \@department\ le notifica a usted que la Comisión 
Revisora de su documento de \@doctype , integrada por los profesores 
\@advisortitle\ \@advisorname\ y \@gradoasesor2\ \@nombreasesor2, ha concluido 
la revisión del mismo y ha dado la aprobación para su respectiva presentación.\\

Por lo anterior, se le autoriza a usted la impresión definitiva de su 
documento de Tesis para la respectiva defensa en el Examen Profesional, a presentarse 
el \@date 
%\vspace{18pt}
%\@author \\

\vspace{60pt}

\begin{center}
Atentamente \\
Zacatecas, Zac., \@fechaaut \\
%\vspace*{12pt}
\vspace*{28pt}                       
\underline{\makebox[2.8in][r]{}}\\                  
\@dirtitle \hspace{3pt} \@dirname \\
Director de la \@department
\end{center}
}

\def\endautorizacion{\vfill \hfill\par\newpage}

%=============================================================================
% HOJA DE APROBACION DE EXAMEN PROFESIONAL 
%=============================================================================

\def\prestitle#1{\gdef\@prestitle{#1}}
\def\presname#1{\gdef\@presname{#1}}
\def\vocalatitle#1{\gdef\@vocalatitle{#1}}
\def\vocalaname#1{\gdef\@vocalaname{#1}}
\def\vocalbtitle#1{\gdef\@vocalbtitle{#1}}
\def\vocalbname#1{\gdef\@vocalbname{#1}}
\def\vocalctitle#1{\gdef\@vocalctitle{#1}}
\def\vocalcname#1{\gdef\@vocalcname{#1}}
\def\vocaldtitle#1{\gdef\@vocaldtitle{#1}}
\def\vocaldname#1{\gdef\@vocaldname{#1}}
\gdef\@prestitle{empty}
\gdef\@presname{empty}
\gdef\@vocalatitle{empty}
\gdef\@vocalaname{empty}
\gdef\@vocalbtitle{empty}
\gdef\@vocalbname{empty}
\gdef\@vocalctitle{empty}
\gdef\@vocalcname{empty}
\gdef\@vocaldtitle{empty}
\gdef\@vocaldname{empty}

\def\aprobacion{
             \thispagestyle{empty}
%                  \addtocounter{page}{-1}
\begin{center}
{\bf\expandafter\uppercase\expandafter{Aprobación de 
Examen Profesional}}\\
\vspace*{12pt}
      \begin{figure}[h]
      \centering
      \includegraphics[width=1.5in]{escudo.eps}\\
      \end{figure}
      \vspace{12pt}
\end{center}
      
Se aprueba por unanimidad el Examen Profesional de \@author\ 
presentado el \@date\ para obtener el título de:
                 
\vspace{12pt}
\begin{center}
      {\bf\expandafter\uppercase\expandafter{\@degree}} \\
                  \vspace{48pt}
                  Jurado: \\
                  \vspace{9pt} 
            %Presidente \\           
                \begin{flushright}
                Presidente:~\@advisortitle\ \@advisorname \vspace{4pt} \underline{\makebox[2.2in][r]{}} ~~~~~~\\
                \end{flushright} 
                  \vspace{9pt}
            %Primer vocal \\                  
                \begin{flushright}                
               Primer vocal:~\@gradoasesor2\ \@nombreasesor2 \vspace{4pt} \underline{\makebox[2.2in][r]{}} ~~~~~~\\
                \end{flushright}
                  \vspace{9pt}         
            %Segundo vocal \\   
                \begin{flushright}                          
                 Segundo vocal:~\@vocalbtitle\ \@vocalbname \vspace{4pt} \underline{\makebox[2.2in][r]{}} ~~~~~~\\                  
                \end{flushright}
                  \vspace{9pt}
            %Tercer vocal \\               
                \begin{flushright}                
                Tercer vocal:~\@vocalctitle\ \@vocalcname \vspace{4pt} \underline{\makebox[2.2in][r]{}} ~~~~~~\\                  
                \end{flushright}                  
                  \vspace{9pt}
            %Cuarto vocal \\              
            \begin{flushright}
             Cuarto vocal:~\@vocaldtitle\ \@vocaldname \vspace{4pt} \underline{\makebox[2.2in][r]{}} ~~~~~~\\                                                           
            \end{flushright}
                \end{center}
           }

\def\endaprobacion{\vfill \hfill\par\newpage}

%=============================================================================
% RESUMEN (En español)
% The ESP abstract should begin with the author and title in a standard format.
% After the author comes the advisor and university. After these lines comes
% a bunch of double spaced text to make up the standard abstract.
% After the abstract, the advisor's approval signature follows.
% This page is not numbered and is delivered seperately to the thesis office.
%=============================================================================

\def\advisortitle#1{\gdef\@advisortitle{#1}}
\def\advisorname#1{\gdef\@advisorname{#1}}
\gdef\@advisortitle{empty}
\gdef\@advisorname{empty}

\def\resumenespa{
\addcontentsline{toc}{chapter}{Resumen}
             \thispagestyle{empty}
%                  \addtocounter{page}{-1}
                \begin{center}
                  %{\bf\expandafter\uppercase\expandafter{\@title}}\\
                  %\vspace{12pt}
                  %\@author \\
                  %\vspace{18pt}
                  %Directores de tesis: \@advisortitle\ \@advisorname\ y
                  %\@gradoasesor2\ \@nombreasesor2 \\
                  %\vspace{24pt}
                  \expandafter\uppercase\expandafter{\bf{Resumen}}
                  \vspace{12pt}
                \end{center}
}

\def\endresumenespa{\vfill \hfill\par\newpage}

%=============================================================================
% ABSTRACT (Resumen en inglés)
% The ENG abstract should begin with the author and title in a standard format.
% After the author comes the advisor and university. After these lines comes
% a bunch of double spaced text to make up the standard abstract.
% After the abstract, the advisor's approval signature follows.
% This page is not numbered and is delivered seperately to the thesis office.
%=============================================================================

%\def\advisortitle#1{\gdef\@advisortitle{#1}}
%\def\advisorname#1{\gdef\@advisorname{#1}}
%\gdef\@advisortitle{empty}
%\gdef\@advisorname{empty}
\def\englishtitle#1{\gdef\@englishtitle{#1}}
\gdef\@englishtitle{empty}

\def\resumeningles{
\addcontentsline{toc}{chapter}{Abstract}
             \thispagestyle{empty}
%                  \addtocounter{page}{-1}
                \begin{center}
                  %{\bf\expandafter\uppercase\expandafter{\@englishtitle}}\\
                  %\vspace{12pt}
                  %\@author \\
                  %\vspace{15pt}
                  %Thesis supervisor: \@advisortitle\ \@advisorname \\
                  %\vspace{25pt}
                  \expandafter\uppercase\expandafter{\bf{Abstract}}
                  \vspace{12pt}
                \end{center}
}
\def\endresumeningl{\vfill \hfill\par\newpage}

%=============================================================================
% RESUMEN
% The abstract should begin with two single-spaced lines describing
% the author and title in a standard format.  After these lines comes
% the standard abstract.
%=============================================================================
%\def\resumen{
%  \chapter*{RESUMEN}
%  \addcontentsline{toc}{chapter}{RESUMEN}
%  \relax\markboth{RESUMEN}{RESUMEN}}
%\def\endresumen{\par\newpage}
%
%=============================================================================
% ABSTRACT
% El abstract debe comenzar con dos renglones a espacio sencillo describiendo
% al título, y autor.  After these lines comes
% the standard abstract.
%=============================================================================
%\def\abstract{
%  \chapter*{Abstract}
%  \addcontentsline{toc}{chapter}{Abstract}
%  \relax\markboth{Abstract}{Abstract}}
%\def\endabstract{\par\newpage}
%
%=============================================================================
% PAGE STYLES
%=============================================================================
% The page style 'foo' is defined by defining the command \ps@foo.  This
% command should make only local definitions.  There should be no stray
% spaces in the definition, since they could lead to mysterious extra
% spaces in the output.
%
% The \ps@... command defines the macros \@oddhead, \@oddfoot,
% \@evenhead, and \@evenfoot to define the running heads and
% feet---e.g., \@oddhead is the macro to produce the contents of the
% heading box for odd-numbered pages.  It is called inside an \hbox of
% width \textwidth.
%
% To make headings determined by the sectioning commands, the page style
% defines the commands \chaptermark, mark, ... , where
% \chaptermark{TEXT} is called by \chapter to set a mark, and so on.
% The \...mark commands and the \...head macros are defined with the
% help of the following macros.  (All the \...mark commands should be
% initialized to no-ops.)
%
% MARKING CONVENTIONS:
% LaTeX extends TeX's \mark facility by producing two kinds of marks
% a 'left' and a 'right' mark, using the following commands:
%     \markboth{LEFT}{RIGHT} : Adds both marks.
%     \markright{RIGHT}      : Adds a 'right' mark.
%     \leftmark  : Used in the \@oddhead, \@oddfoot, \@evenhead or \@evenfoot
%                  macro, gets the current 'left'  mark.  Works like TeX's
%                  \botmark command.
%     \rightmark : Used in the \@oddhead, \@oddfoot, \@evenhead or \@evenfoot
%                  macro, gets the current 'right'  mark.  Works like TeX's
%                  \firstmark command.
% The marking commands work reasonably well for right marks 'numbered
% within' left marks--e.g., the left mark is changed by a \chapter command and
% the right mark is changed by a \section command.  However, it does
% produce somewhat anomalous results if two \bothmark's occur on the same page.
%=============================================================================
\mark{{}{}}   % Initializes TeX's marks

%=============================================================================
% Definition of 'headings' page style
%  Note use of ##1 for parameter of \def\chaptermark inside \def\ps@headings.
%=============================================================================
\if@twoside                                                  % If two-sided
  \def\ps@headings{
    \def\@oddfoot{}\def\@evenfoot{}         % No feet.
    \def\@evenhead{\rm \thepage\hfil \sl \leftmark}          % Left heading.
    \def\@oddhead{\hbox{}\sl \rightmark \hfil \rm\thepage}   % Right heading.
    \def\chaptermark##1{\markboth {\uppercase{\ifnum \c@secnumdepth >\m@ne
      \@chapapp\ \thechapter. \ \fi ##1}}{}}%
    \def\sectionmark##1{\markright {\uppercase{\ifnum \c@secnumdepth >\z@
      \thesection. \ \fi ##1}}}}
\else                                                        % If one-sided
  \def\ps@headings{
    \def\@oddfoot{}\def\@evenfoot{}                          % No feet.
    \def\@oddhead{\hbox {}\sl \rightmark \hfil \rm\thepage}  % Heading.
    \def\chaptermark##1{\markright {\uppercase{\ifnum \c@secnumdepth >\m@ne
      \@chapapp\ \thechapter. \ \fi ##1}}}}
\fi

%=============================================================================
% Definition of 'thesis' page style -- riedl
%  -- no footer
%  -- header is just page number, in upper right corner
%  -- chaptermark and sectionmark are not needed
%  Note use of ##1 for parameter of \def\chaptermark inside \def\ps@headings.
%=============================================================================
\if@twoside                                                  % If two-sided
  \def\ps@thesis{
    \def\@oddfoot{}\def\@evenfoot{}                          % No feet.
    \def\@evenhead{\rm \thepage\hfil}                        % Left heading.
    \def\@oddhead{\hbox{}\hfil\rm\thepage}                   % Right heading.
    \def\chaptermark##1{}
    \def\sectionmark##1{}}
\else                                                        % If one-sided
  \def\ps@thesis{
    \def\@oddfoot{}\def\@evenfoot{}                          % No feet.
    \def\@oddhead{\hbox{}\hfil\rm\thepage}                   % Heading.
    \def\chaptermark##1{}}
\fi

%=============================================================================
% The following code sets \@currentipfile to current include file.
% taken from chapterbib.sty by
% Niel Kempson (JANET: rmcs-tex@uk.ac.cranfield.cdvc)

\let\@debugoutput=\@gobble

%-------------------------------------------------------------------------
% Keep track of files currently opened.
% 1. The counter @ipfilectr keeps track of the currently open I/P file,
%    and @includefilectr keeps track of the \include files.
% 2. @ipfilectr is set to @includefilectr, except when processing the root
%    file. An arbitrary large (number greater than the maximum number of
%    include files) is used to denote the root file.
%-------------------------------------------------------------------------
\newcounter{@ipfilectr}                 % current I/P file,
\def\the@ipfilectr{\roman{@ipfilectr}}  % displayed as roman numerals
\newcounter{@includefilectr}            % current \include file
\def\@rootfilevalue{1988}               % arbitrary value for root file
\def\savedjobname{\jobname}             % save root file name as a macro

%
%-----------------------------------------------------------------------
% The name of the current I/P file is stored in macro \@currentipfile
% and the global macro ipfile@the@ipfilectr is set to this name whenever
% \@currentipfile changes.
%-----------------------------------------------------------------------
\def\@currentipfile{\jobname}           % initialize to root file
\setcounter{@ipfilectr}{\@rootfilevalue}        % ditto
\global\@namedef{ipfile@\the@ipfilectr}{\@currentipfile}

%
%-----------------------------------------------------------------------
% Redefine the \include macro (taken from LATEX.TEX of 25-JAN-1988) so
% that when a file is \include'd
% 1. \@currentipfile is updated and \includefilectr incremented
% 2. ipfilectr set to includefilectr
% 3. global macro ipfile@\the@ipfilectr set to \@currentfile
%
% When the \include'd file is finished with
% 1. reset \@currentipfile to the root file name (\jobname)
% 2. reset ipfilectr to \@rootfilevalue
% 3. global macro ipfile@\the@ipfilectr set to \@currentfile
%-----------------------------------------------------------------------
\def\include#1{\clearpage
\def\@currentipfile{#1}
\@debugoutput{\string\@currentipfile = `\@currentipfile'}
\addtocounter{@includefilectr}{1}
\setcounter{@ipfilectr}{\value{@includefilectr}}
\global\@namedef{ipfile@\the@ipfilectr}{\@currentipfile}
%
\if@filesw \immediate\write\@mainaux{\string\@input{#1.aux}}\fi
\@tempswatrue\if@partsw \@tempswafalse\def\@tempb{#1}\@for
\@tempa:=\@partlist\do{\ifx\@tempa\@tempb\@tempswatrue\fi}\fi
\if@tempswa \if@filesw \let\@auxout=\@partaux
\immediate\openout\@partaux #1.aux
\immediate\write\@partaux{\relax}\fi\@input{#1.tex}\clearpage
\@writeckpt{#1}\if@filesw\immediate\closeout\@partaux \fi
\let\@auxout=\@mainaux\else\@nameuse{cp@#1}\fi
%
\def\@currentipfile{\jobname}
\@debugoutput{\string\@currentipfile = `\@currentipfile'}
\setcounter{@ipfilectr}{\@rootfilevalue}
\global\@namedef{ipfile@\the@ipfilectr}{\@currentipfile}}
%=============================================================================

%=============================================================================
% Definition of 'thesisdraft' page style -- riedl --dfk
%  -- header is just page number, in upper right corner
%  -- chaptermark and sectionmark are not needed
%  -- footer is set to identify this as a draft version (timestamped)
%  -- footer may be changed by changing \@draftfoot
%  Note use of ##1 for parameter of \def\chaptermark inside \def\ps@headings.
%=============================================================================
\newcounter{currenthour} \newcounter{currentminute}
\newcommand{\now}{
  \setcounter{currenthour}{\time}\divide \c@currenthour by 60
  \setcounter{currentminute}{\time}
    \multiply \c@currenthour by 60
    \addtocounter{currentminute}{-\c@currenthour}
    \divide \c@currenthour by 60
  \ifnum \c@currenthour=0 0 \fi    % force at least one digit in hour
  \thecurrenthour:%
  \ifnum \c@currentminute<10 0\fi % force two digits in minute
  \thecurrentminute
}

% changed by mccauley to print name of currently included file.
\def\@draftfooter{
  {\rm DRAFT: Do Not Distribute}\hfil{\now\ \today}
    \hfil{\tt \@currentipfile.tex}}

\if@twoside                                                  % If two-sided
  \def\ps@thesisdraft{
    \def\@oddfoot {\@draftfooter}                            % Footer
    \def\@evenfoot{\@draftfooter}
    \def\@evenhead{\rm \thepage\hfil}                        % Left heading.
    \def\@oddhead{\hbox{}\hfil\rm\thepage}                   % Right heading.
    \def\chaptermark##1{}
    \def\sectionmark##1{}}
\else                                                        % If one-sided
  \def\ps@thesisdraft{
    \def\@oddfoot {\@draftfooter}                            % Footer
    \def\@evenfoot{}
    \def\@oddhead{\hbox{}\hfil\rm\thepage}                   % Heading.
    \def\chaptermark##1{}}
\fi

%=============================================================================
% Definition of 'myheadings' page style.
%=============================================================================
\def\ps@myheadings{
  \def\@oddhead{\hbox{}\sl\rightmark \hfil \rm\thepage}
  \def\@oddfoot{}
  \def\@evenhead{\rm \thepage\hfil\sl\leftmark\hbox {}}
  \def\@evenfoot{}
  \def\sectionmark##1{}
  \def\subsectionmark##1{}}

%=============================================================================
% MISCELANEA
%=============================================================================
%\def\mybox{\fbox{\hbox to 429bp {\vbox to 630bp {\null}}}}

%\def\pagebox{ \setlength{\unitlength}{1.0bp}
%             \begin{picture}(0,0)(36,558) \put(0,0){\mybox} \end{picture}
%              \begin{picture}(0,0)(18,523.25) \put(0,0){\mybox} \end{picture}
%              \begin{picture}(0,0)(18,500.5) \put(0,0){\mybox} \end{picture}
%           }


\def\draftmargins{
  \if@msmargin % Draw box for 1.25" left margins
      \special{!userdict begin /eop-hook{
            gsave 0.5 setlinewidth newpath
            90 72 moveto 540 72 lineto 540 720 lineto 90 720 lineto
            90 72 lineto
            stroke grestore
%        Cross at 0.5" down and 1" from right edge
%           gsave 0.5 setlinewidth newpath
%           522 756 moveto 558 756 lineto
%           540 738 moveto 540 774 lineto
            stroke grestore
           }def end
   }
  \else % Draw box for 1" margins
      \special{!userdict begin /eop-hook{
            gsave 0.5 setlinewidth newpath
            72 72 moveto 540 72 lineto 540 720 lineto 72 720 lineto
            72 72 lineto
            stroke grestore
%        Cross at 0.5" down and 1" from right edge
%           gsave 0.5 setlinewidth newpath
%           522 756 moveto 558 756 lineto
%           540 738 moveto 540 774 lineto
            stroke grestore
           }def end
   }

  \fi
}

\def\draftscreen{
   \special{!userdict begin /bop-hook{gsave
            260 185 translate
            65 rotate /Times-Roman findfont 144 scalefont setfont
            0 0 moveto 0.9 setgray (DRAFT) show grestore
            }def end
   }
}

%============================================================================
% VERBATIMFILE
%============================================================================
% \verbatimfile{<filename>}    for verbatim inclusion of a file
% - Note that the precise layout of line breaks in this file is important!
% - added the \singlespace - EB
%============================================================================
\def\verbatimfile#1{\begingroup \singlespace
                    \@verbatim \frenchspacing \@vobeyspaces
                    \input#1 \endgroup
}


%=============================================================================
% SEPARATOR Pages
%   Creates a blank page with a text centered horizontally and vertically.
%   The page is neither counted nor numbered.
%   These pages are required in the thesis format before sections such
%   as appendices, vita, bibliography, etc.
%=============================================================================
\def\separatorpage#1{
  \newpage
  \thispagestyle{empty}
  \addtocounter{page}{-1}
  \null
  \vfil\vfil
  \begin{center}
    {\bface #1}
  \end{center}
  \vfil\vfil
  \newpage}

%=============================================================================
% \@chapapp is initially defined to be 'Chapter'.  The \appendix
% command redefines it to be 'Appendix'.
%=============================================================================
\def\@chapapp{Capítulo}

%=============================================================================
% APÉNDICES
%=============================================================================
% The \appendix command must do the following:
%    -- reset the chapter counter to zero
%    -- set \@chapapp to Appendix (for messages)
%    -- redefine the chapter counter to produce appendix numbers
%    -- reset the section counter to zero
%    -- redefine the \chapter command if appendix titles and headings are
%       to look different from chapter titles and headings.
%    -- produce special toc, lot and lof entries
%=============================================================================
\def\appendix{
  \@appendix{Apéndice}
  \def\@chapstyle{appendix}
  \def\@chapnumber{Apéndice}
  \let\@chapterhead\@appendixchapterhead
  \let\@schapterhead\@appendixchapterhead}
\def\endappendix{
  \def\@chaptername##1{\uppercase {##1}}
  \def\@chapstyle{chapter}
  \def\@chapnumber{\thechapter.}
  \let\@chapterhead\@makechapterhead
  \let\@schapterhead\@makeschapterhead}

\def\l@appendix#1#2{\pagebreak[3]
  \vskip\@tocchapterskip % space above chapter line
  \@dottedtocline{0}{0pt}{7.3em}{{\bface #1}}{#2}}

\def\appendices{
  \addtocontents{toc}{
    \protect\addvspace{\@tocchapterskip}
    \hbox to \hsize{\bface Apéndices \hfil}
    \protect\addvspace{\@tocchapterskip}}
  \@appendix{Apéndices}
  \def\@chaptername##1{##1}
  \def\@chapstyle{appendices}
  \def\@chapnumber{Apéndice\ \thechapter:}
  \let\@chapterhead\@appendiceschapterhead
  \let\@schapterhead\@appendiceschapterhead
  \if@bibpresent\relax \else \separatorpage{}\fi % This
  % is a dirty fix.  For some reason, the TOC isn't
  % properly written unless the simple chapter with
  % sections is written.  The bibliography fixes this
  % problem when its present.  If its not, then the
  % seperator page is inserted to fix things.  EB
}

\let\endappendices\endappendix

\def\l@appendices#1#2{\@dottedtocline{0}{1.5em}{6.5em}{#1}{#2}}

\def\noappendixtables{\def\@lotheading{}}
\def\noappendixfigures{\def\@lofheading{}}
\def\@lotheading{
  \addtocontents{lot}{
    \protect\stepcounter{listofheading}
    \protect\addvspace{\@tocchapterskip}
    \hbox to \hsize{Apéndice \hfil}
    \hbox to \hsize{Tabla \hfil}}}
\def\@lofheading{
  \addtocontents{lof}{
    \protect\stepcounter{listofheading}
    \protect\addvspace{\@tocchapterskip}
    \hbox to \hsize{Apéndice \hfil}
    \hbox to \hsize{Figura \hfil}}}

\def\@appendix#1{\par
  \addtocontents{toc}{{\protect\setcounter{tocdepth}{0}}}
  \setcounter{chapter}{0}
  \setcounter{section}{0}
  \def\@chapapp{APPENDIX}
  \def\thechapter{\Alph{chapter}}
  \@lotheading
  \@lofheading
  %\separatorpage{\uppercase{#1}} 
  % UW doesn't require a seperator page - EB
}

%=============================================================================
% COPYRIGHTPAGE
%=============================================================================
% El ambiente copyright hace lo siguiente:
% - inicia una nueva página sin número de página
% - coloca el texto de copyright ten la parte de abajo.
%=============================================================================
\def\copyrightpage{
  \newpage
  \thispagestyle{empty}      % Sin número de página
  \addtocounter{page}{-1}
  \chapter*{}                % Requiredo para que trabaje \vfill
  \begin{center}
   \vfill
   \copyright\ Copyright de \@author\ \@date\\
   Derechos reservados
  \end{center}}

%=============================================================================
% GLOSARIO
%=============================================================================
% El ambiente glosario hace lo siguiente:
% - genera espacio en el índice de Contenidos para Glosario
% - comienza una nueva página con el encabezado Glosario
%   (2" abajo del márgen superior) y el texto correspondiente
%=============================================================================
\def\glosario{
  \chapter*{Glosario}
  \addcontentsline{toc}{chapter}{Glosario}}
\def\endglosario{\par\newpage}

%=============================================================================
% NOMENCLATURA
%=============================================================================
% El ambiente nomenclatura hace lo siguiente:
% - genera espacio en el índice de contenidos para Nomenclatura
% - comienza una nueva página con el encabezado Nomenclatura
%   (2" abajo del márgen superior) y el texto correspondiente
%=============================================================================
\def\nomenclatura{ %\separatorpage{DISCARD THIS PAGE}
  \chapter*{Nomenclatura}
  \addcontentsline{toc}{chapter}{Nomenclatura}}
\def\endnomenclatura{\par\newpage}

%=============================================================================
% LISTA DE SIMBOLOS
%=============================================================================
% El ambiente listofsymbols hace lo siguiente:
% - genera espacio en el índice de contenidos para Lista de símbolos
% - comienza una nueva página con el encabezado Lista de símbolos
%   (2" abajo del márgen superior) y el texto correspondiente
%=============================================================================
\def\listofsymbols{ %\separatorpage{DISCARD THIS PAGE}
  \chapter*{Lista de símbolos}
  \addcontentsline{toc}{chapter}{Lista de símbolos}}
\def\endlistofsymbols{\par\newpage}

%=============================================================================
% CURRICULUM VITAE
%=============================================================================
% El ambiente vita hace lo siguiente:
% - genera espacio en el índice de Contenidos para Curriculum Vitae
% - comienza una nueva página con el encabezado Curriculum Vitae
%   (2" abajo del márgen superior) y el texto correspondiente
%=============================================================================
\def\vita{
%  \separatorpage{Curriculum Vitae}         % UW doesn't require this EB
  \chapter*{Curriculum Vitae}
  \addcontentsline{toc}{chapter}{Curriculum Vitae}}
\def\endvita{\par\newpage}

%=============================================================================
% DEDICATORIA
%=============================================================================
% El ambiente dedicatoria hace lo siguiente: 
% - inicia una nueva página
% - centra el texto verticalmente
% - coloca el texto en ambiente centrado
%=============================================================================
\def\dedicatoria{
  \newpage
  \null\vfil
  \begin{center}}
\def\enddedicatoria{\end{center}\par\vfil\newpage}

%=============================================================================
% AGRADECIMIENTOS
%=============================================================================
% El ambiente agradecimientos:
% - inicia una nueva página con el encabezado Agradecimientos 
%   (2" abajo del márgen superior) y el texto correspondiente 
%   en forma de capítulo 
%=============================================================================
\def\agradecimientos{
  \chapter*{Agradecimientos}
  \doublespace}
\def\endaagradecimientos{\par\newpage}

%=============================================================================
% FECHAS
%=============================================================================
%\def\today{\ifcase\month\or
  %January\or February\or March\or April\or May\or June\or
  %July\or August\or September\or October\or November\or December\fi
  %\space\number\day, \number\year}
\newcount\@testday
\def\today{\@testday=\day
  \ifnum\@testday>30 \advance\@testday by -30
  \else\ifnum\@testday>20 \advance\@testday by -20
  \fi\fi
  \number\day\ \
  \ifcase\month\or
    January \or February \or March \or April \or May \or June \or
    July \or August \or September \or October \or November \or December
    \fi\ \number\year
}

%  Single counter for theorems and theorem-like environments:
\newtheorem{theorem}{Theorem}[chapter]
\newtheorem{assertion}[theorem]{Assertion}
\newtheorem{claim}[theorem]{Claim}
\newtheorem{conjecture}[theorem]{Conjecture}
\newtheorem{corollary}[theorem]{Corollary}
\newtheorem{definition}[theorem]{Definition}
\newtheorem{example}[theorem]{Example}
\newtheorem{figger}[theorem]{Figure}
\newtheorem{lemma}[theorem]{Lemma}
\newtheorem{prop}[theorem]{Proposition}
\newtheorem{remark}[theorem]{Remark}

%=============================================================================
% EQUATION and EQNARRAY -- put here because it must follow \chapter definition
%
% \newcounter{equation}
%=============================================================================
% \jot = 3pt      % Extra space added between lines of an eqnarray environment
\@addtoreset{equation}{chapter}   % Makes \chapter reset 'equation' counter.
\def\theequation{{\rm \thechapter}.\arabic{equation}}

%=============================================================================
% The macro \@eqnnum defines how equation numbers are to appear in equations.
%=============================================================================
\def\@eqnnum{(\theequation)}

%=============================================================================
% INITIALIZATION   -  Default initializations
%=============================================================================
\ps@thesis                  % Estilo 'plain' de página
\pagenumbering{arabic}      % Números arábicos de página
\onecolumn                  % Columna sencilla.
\raggedbottom               % Ragged bottom
%\if@twoside\else\raggedbottom\fi % Ragged bottom unless twoside option.

\global\newif\if@bibpresent   % Flag for the appendicies fix.  If bib present
                              % then the seperator page isn't needed...

%-----------------------------------------------------------------------------
% Define \bface to be \bf so that all default boldface text can be
% disabled by redefining \bface to be \null.  -dfk
%-----------------------------------------------------------------------------
\newcommand{\bface}{\bf}

%-----------------------------------------------------------------------------
% The following penalties are for conforming with Purdue's
% requirements for the thesis format.
%-----------------------------------------------------------------------------
\brokenpenalty=10000
\clubpenalty=10000
\widowpenalty=10000

%=============================================================================
% FIN DEL ARCHIVO
%=============================================================================
