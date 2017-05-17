% =========================================================================
% Class for the smartplus documentation
% 

\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{clsbook}[2017/05/16 smartplus documentation book class]

%% Les differents styles 
\newif\if@stylebook \@stylebookfalse

%% class book, with various options
\PassOptionsToClass{USletter}{book}
\PassOptionsToClass{11pt}{book}
\PassOptionsToClass{twoside}{book}
\DeclareOption*{\PassOptionsToClass{\CurrentOption}{book}}
\DeclareOption{stylebook}{\@stylebooktrue}
\DeclareOption*{}
\ProcessOptions
\LoadClass{book}

% ==========================================================================
% PACKAGES STANDARD REQUIS

%fonts and style
\RequirePackage[latin1]{inputenc}
\RequirePackage[T1]{fontenc}
\RequirePackage[francais]{babel}
\frenchspacing
\RequirePackage{eurosym}
\RequirePackage{xcolor}
\RequirePackage{gensymb}
\RequirePackage{lmodern}

%selection of the stylebook
\if@stylebook
  \RequirePackage{stylebook}
\fi
