;;-----------------------------------------------------------+
;;                                                           |
;; EMACS APOLLO STYLES MODE                                  |
;;                                                           |
;;-----------------------------------------------------------+
;;                                                           |
;;  AUTHOR: Jamie Estill                                     |
;; CONTACT: jestill_@_sourceforge.org                        |
;; STARTED: 12/12/2006                                       |
;; UPDATED: 12/12/2006                                       |
;;   NOTES:                                                  |
;;    Modified from the apollo tiers mode file.              |
;;                                                           |
;;-----------------------------------------------------------+

(defvar apollo-styles-mode-hook nil)

;;-----------------------------------------------------------+
;; FILE EXTENSIONS THIS APPLIES TO                           |
;;-----------------------------------------------------------+
(add-to-list 'auto-mode-alist '("\\.style\\'" . apollo-styles-mode))

;;-----------------------------------------------------------+
;; FONT-LOCK-KEYWORDS                                        |
;; 11/10/2006                                                |
;;-----------------------------------------------------------+
(defconst apollo-styles-font-lock-keywords
  (list
   ;; The regexp-opt function optimizes the regular expression 
   ;; for a list of characters 

   ;;-----------------------------+
   ;; VARIABLE NAMES              |
   ;;-----------------------------+
   ;; In Apollo these are refered to as keys
   (cons (concat "\\<\\(" (regexp-opt '
			   ("ImportStyle" 
			    "Types" 
			    "Chromosomes"
			    "Draw3D" 
			    "ShowAnnotations" 
			    "ShowResults"
			    "EnableEditing" 
			    "EnableNavigationManager"
			    "AnnotationBackgroundColor"
			    "FeatureBackgroundColor" 
			    "EdgematchColor"
			    "EdgematchWidth" 
			    "FeatureLabelColor"
			    "AnnotationLabelColor" 
			    "CoordBackgroundColor"
			    "CoordForegroundColor" 
			    "UserTranscriptColouring"
			    "SelectionColor" 
			    "DashSets"
			    "SequenceGapColor" 
			    "DrawOutline"
			    "DefaultFeatureLabelFont" 
			    "DefaultFont"
			    "ExonDetailEditorSequenceFont"
			    "ExonDetailEditorBackgroundColor1"
			    "ExonDetailEditorBackgroundColor2"
			    "ExonDetailEditorFeatureColor1"
			    "ExonDetailEditorFeatureColor2"
			    "NameAdapterInstall" 
			    "ShowIdField"
			    "TranscriptSymbolEditable"
			    "ShowIsDicistronicAnnotInfoCheckbox"
			    "ShowEvaluationOfPeptide"
			    "ShowIsProblematicAnnotInfoCheckbox"
			    "ShowFinishedAnnotInfoCheckbox"
			    "ShowDBCrossRefAnnotInfoTable"
			    "ShowReplaceStopAnnotInfoCheckbox"
			    "EnableTranslationalFrameShiftEditing"
			    "EnableSequenceErrorEditing"
			    "ShowOwnershipAnnotMenuItem"
			    "ShowTranscriptFinishedAnnotMenuItem"
			    "DisplayPreferences"
			    ;; URLS
			    "ExternalRefURL" 
			    "DatabaseURLField"
			    "GeneURL" 
			    "BandURL" 
			    "ScaffoldURL"
			    "RangeURL"
			    "SequenceURL"
			    ;; May be deprecated
			    "OverlapDefinition"
			    ;; Other Variables
			    "UserInfo"
			    "PeptideStatus"
			    "ResultTag"
			    "EnableIgbHttpConnection" 
			    "DatabaseList"
			    "AnnotationComment" 
			    "TranscriptComment") 
			   t) "\\)\\>")'font-lock-variable-name-face)


   )
"Minimal highlighting expression for apollo tiers mode")

;;-----------------------------------------------------------+
;; COMMENTS IN APOLLO STYLES FILES                           |
;;-----------------------------------------------------------+
(defvar apollo-styles-mode-syntax-table
  (let (( apollo-styles-mode-syntax-table (make-syntax-table)))

    ;; < is the open comment
    ;; > is the close comment
    ;; The following thinks that / is a comment
    ;;(modify-syntax-entry ?/ "<" apollo-styles-mode-syntax-table)
    (modify-syntax-entry ?/ ". 12" apollo-styles-mode-syntax-table)
    (modify-syntax-entry ?\n ">" apollo-styles-mode-syntax-table)
    apollo-styles-mode-syntax-table)
  "Syntax table for apollo-tiers-mode")

;;-----------------------------------------------------------+
;; APOLLO TIERS MODE FUNCTION                                |
;;-----------------------------------------------------------+
;; This is the function that will be called by Emacs when 
;; the mode is started
(defun apollo-styles-mode ()
  "Major mode for editing apollo tiers files"
  (interactive)
  (kill-all-local-variables)
  (set-syntax-table apollo-styles-mode-syntax-table)
  (set (make-local-variable 'font-lock-defaults) 
       '(apollo-styles-font-lock-keywords))
  (setq mode-name "Apollo Style")
  (run-hooks 'apollo-styles-mode-hook)
)

(provide 'apollo-styles-mode)