;;-----------------------------------------------------------+
;;                                                           |
;; EMACS APOLLO TIERS MODE                                   |
;;                                                           |
;;-----------------------------------------------------------+
;;                                                           |
;;  AUTHOR: Jamie Estill                                     |
;; CONTACT: jestill_@_sourceforge.org                        |
;; STARTED: 5/26/2006                                        |
;; UPDATED: 11/10/2006                                       |
;;   NOTES:                                                  |
;;    Comments and URLS fight over precedence for syntax     |
;;    highlighting.                                          |
;;                                                           |
;;-----------------------------------------------------------+

(defvar apollo-tiers-mode-hook nil)

;;-----------------------------------------------------------+
;; FILE EXTENSIONS THIS APPLIES TO                           |
;;-----------------------------------------------------------+
(add-to-list 'auto-mode-alist '("\\.tiers\\'" . apollo-tiers-mode))

;;-----------------------------------------------------------+
;; FONT-LOCK-KEYWORDS                                        |
;; 11/10/2006                                                |
;;-----------------------------------------------------------+
(defconst apollo-tiers-font-lock-keywords
  (list
   ;; The regexp-opt function optimizes the regular expression 
   ;; for a list of characters 
   ;;-----------------------------+
   ;; COMMENTS                    |
   ;;-----------------------------+
   ;; The following works okay, but the use of // as comments takes
   ;; precedence even when using fon-lock-comment-face
   ;;'("\/\/.*\n" . font-lock-comment-face)
   ;; Changed the above to the following to only match at the
   ;; beginning of a line.
   '("^\/\/.*\n" . font-lock-comment-face)
   ;; The following comments out lines starting with #
   ;;'("\#.*\n" . font-lock-comment-face)
   ;; Changed from the above to the following
   '("^\#.*\n" . font-lock-comment-face)

   ;;-----------------------------+
   ;; VARIABLE NAMES              |
   ;;-----------------------------+
   (cons (concat "\\<\\(" (regexp-opt '
			   ("tiername" "datatype" "glyph" "utr_color" 
			    "name_method" "overlap_method" "idformat" "color" 
			    "label" "curated" "groupby" "number_of_levels" 
			    "warnonedit" "sorted" "expanded" "usescore" 
			    "freshdate" "visible" "minscore" "maxscore"
			    "annot_type" "weburl" "labeled" "sortbycolumn" 
			    "column" "maxrows" "reversesort"
			    "typename" "resulttype" "scorethreshold") 
			   t) "\\)\\>")'font-lock-variable-name-face)


   ;;-----------------------------+
   ;; URLS                        |
   ;;-----------------------------+
   ;; The following works okay, but the use of // as comments takes
   ;; precedence even when using fon-lock-comment-face
   '("\\http.*\n" . font-lock-keyword-face)


   ;;-----------------------------+
   ;; COLUMNS THAT WILL BE SHOWN  |
   ;;-----------------------------+
   (cons (concat "\\<\\(" (regexp-opt '
			   ("GENOMIC_RANGE" "ID" "GENOMIC_LENGTH" 
			    "MATCH_RANGE" "MATCH_LENGTH" "SCORE" "query_frame" 
			    "HOMOLOGY" "NAME" "identity" "expect" 
			    "probability" "score") 
			   t) "\\)\\>")'font-lock-string-face)
   ;;-----------------------------+
   ;; BOOLEANS                    |
   ;;-----------------------------+
   (cons (concat "\\<\\(" (regexp-opt '
			   ("true" "false")
			   t) "\\)\\>")'font-lock-keyword-face)
   ;;-----------------------------+
   ;; GLYPHS                      |
   ;;-----------------------------+
   (cons (concat "\\<\\(" (regexp-opt '
			   ("PromoterGlyph" "SiteCodon" 
			    "DrawableResultFeatureSet" "DrawableGeneFeatureSet"
			    "DoubleHeadedArrow" "Triangle" "Zigzag" 
			    "ThinRectangle") 
			   t) "\\)\\>")'font-lock-constant-face)
   ;;-----------------------------+
   ;; GROUP-BY METHODS            |
   ;;-----------------------------+
   (cons (concat "\\<\\(" (regexp-opt '
			   ("GENE" "SINGLE" "HOMOLOGY") 
			   t) "\\)\\>")'font-lock-type-face)
   ;;-----------------------------+
   ;; OVERLAP METHODS             |
   ;;-----------------------------+
   (cons (concat "\\<\\(" (regexp-opt '
			   ("SimpleOverlap" "ORF_Overlap" "NoOverlap") 
			   t) "\\)\\>")'font-lock-type-face)
   ;;-----------------------------+
   ;; THE NAME ADAPTERS           |
   ;;-----------------------------+
   (cons (concat "\\<\\(" (regexp-opt '
			   ("SimpleNameAdapter" "FlyNameAdapter" 
			    "DefaultNameAdapter" "RiceNameAdapter"
			    "GmodNameAdapter" "TigrSybilNameAdapter") 
			   t) "\\)\\>")'font-lock-function-name-face)

   ;;-----------------------------+
   ;; MAIN HEADERS                |
   ;;-----------------------------+
   '("\\[Tier\\]" . font-lock-warning-face)
   '("\\[Type\\]" . font-lock-warning-face)


   )
"Minimal highlighting expression for apollo tiers mode")

;;-----------------------------------------------------------+
;; COMMENTS IN APOLLO TIERS FILES                            |
;;-----------------------------------------------------------+
;; Commented out the following - 11/10/2006
;; This should only show commenting for lines that start
;; out as comments
;;(defvar apollo-tiers-mode-syntax-table
;;  (let (( apollo-tiers-mode-syntax-table (make-syntax-table)))
;;    (modify-syntax-entry ?/ ". 124b" apollo-tiers-mode-syntax-table)
;;    (modify-syntax-entry ?\n "> b" apollo-tiers-mode-syntax-table)
;;    ;; The following two appear to be borked
;;    ;;(modify-syntax-entry ?# "<"  apollo-tiers-mode-syntax-table)
;;    ;;(modify-syntax-entry ?\n ">"  apollo-tiers-mode-syntax-table)
;;    apollo-tiers-mode-syntax-table)
;;"Syntax table for apollo-tiers-mode")


;;-----------------------------------------------------------+
;; APOLLO TIERS MODE FUNCTION                                |
;;-----------------------------------------------------------+
;; This is the function that will be called by Emacs when 
;; the mode is started
(defun apollo-tiers-mode ()
  "Major mode for editing apollo tiers files"
  (interactive)
  (kill-all-local-variables)
;; Commented out the following 11/10/2006
;;  (set-syntax-table apollo-tiers-mode-syntax-table)
  (set (make-local-variable 'font-lock-defaults) 
       '(apollo-tiers-font-lock-keywords))
  (setq mode-name "Apollo Tiers")
  (run-hooks 'apollo-tiers-mode-hook)
)

(provide 'apollo-tiers-mode)