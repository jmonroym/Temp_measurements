\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{unlsilabsop}[2014/06/27 v0.0 Modified article class for standard operating procedures in UNL SiLab]
% http://stackoverflow.com/questions/581916/how-do-you-extend-article-document-class-in-latex

% Passes and class options to the underlying article class
\DeclareOption*{\PassOptionsToClass{\CurrentOption}{article}}
\ProcessOptions
\LoadClass{article}

\RequirePackage[utf8]{inputenc}
% Redefine the page margins
\RequirePackage[left=1in,right=1in,top=1in,bottom=1in]{geometry}

% Modifications to the section titles
\RequirePackage{titlesec}
\RequirePackage{graphicx}
\RequirePackage{tocloft}% http://ctan.org/pkg/tocloft
\RequirePackage{float}
\restylefloat{table}
\setlength{\cftsecnumwidth}{4em}% Set length of number width in ToC for \subsection
\setlength{\cftsubsecnumwidth}{4em}% Set length of number width in ToC for \subsection

\renewcommand{\thesection}{\Roman{section}}
\titleformat{\section}{\normalfont\bfseries}
  {\makebox[3em][l]{\thesection{}.}}{0pt}{}
  \titleformat{\subsection}{\normalfont\bfseries}
    {}{0pt}{}

    % Modification of title block
    \RequirePackage{titling}
    \RequirePackage{multirow}
    \RequirePackage{booktabs}
    \newcommand{\sopid}[1]{\newcommand{\thesopid}{#1}}
    \newcommand{\sopversion}[1]{\newcommand{\thesopversion}{#1}}
    \newcommand{\sopabstract}[1]{\newcommand{\thesopabstract}{#1}}
    \newcommand{\approved}[1]{\newcommand{\theapproved}{#1}}
    % Ref: http://tex.stackexchange.com/questions/3988/titlesec-versus-titling-mangling-thetitle
    \let\oldtitle\title
    \renewcommand{\title}[1]{\oldtitle{#1}\newcommand{\mythetitle}{#1}}
    \renewcommand{\maketitle}{%
	\hspace*{-0.7cm}\begin{tabular}{|c|p{2.5cm}p{10.0cm}|} \hline
	    \multirow{5}{*}{\includegraphics[width=2.5cm]{CLASSElogo.png}} & \multicolumn{2}{c|}{\centering \sffamily Standard Operating Procedure} \\
        &          & \\
        &          & \\
	& SOP No.: & \thesopid \\
	& Title: & {\bfseries \mythetitle} \\
CMS Pixel	& Revision: & \thesopversion \\
Phase-II upgrade	& Date: & \thedate \\ \hline
	    \end{tabular}

	\bigskip

    \begin{abstract}\thesopabstract\end{abstract}

	\bigskip

    \tableofcontents
    \newpage
    }

% For "Page N of M"
\RequirePackage{lastpage}
% For easier construction of page headers/footers
\RequirePackage{fancyhdr}
\fancypagestyle{plain}{ % for first page
    \fancyhf{}
    \fancyfoot[L]{\framebox{Author: \theauthor}\\ \jobname{}.tex} 
    \fancyfoot[R]{\framebox{Page: \thepage{} of \pageref*{LastPage}}}
    \renewcommand{\headrulewidth}{0pt}
    \renewcommand{\footrulewidth}{0pt}
}

\pagestyle{fancy} % for other pages
\fancyhf{}
\fancyhead[R]{%
    \begin{tabular}{|c|c|c|} \hline %
	SOP No.: & Revision: & Date: \\
	\thesopid & \thesopversion & \thedate \\ \hline
	\end{tabular}%
}
\fancyfoot[L]{\framebox{Author: \theauthor}}
\fancyfoot[R]{\framebox{Page: \thepage{} of \pageref*{LastPage}}} % \pageref* if we use hyperref, \pageref otherwise
\renewcommand{\headrulewidth}{0pt}
\renewcommand{\footrulewidth}{0pt}

% For easier customization of itemized, enumerated, and other lists
\RequirePackage{enumitem}
% For hyperlinked cross-references
\RequirePackage{hyperref}
% Ensure first page is correct style
\thispagestyle{plain}
% That's all, folks!
\endinput
