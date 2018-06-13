;;-----------------------------------------------------------+
;;                                                           |
;; EMACS APOLLO CONFIG FILE MODE                             |
;;                                                           |
;;-----------------------------------------------------------+
;;                                                           |
;;  AUTHOR: Jamie Estill                                     |
;; CONTACT: jestill_@_sourceforge.org                        |
;; STARTED: 12/13/2006                                       |
;; UPDATED: 12/13/2006                                       |
;;   NOTES:                                                  |
;;    Modified from the apollo tiers mode file.              |
;;                                                           |
;;-----------------------------------------------------------+

(defvar apollo-config-mode-hook nil)

;;-----------------------------------------------------------+
;; FILE EXTENSIONS THIS APPLIES TO                           |
;;-----------------------------------------------------------+
(add-to-list 'auto-mode-alist '("\\.cfg\\'" . apollo-config-mode))

;;-----------------------------------------------------------+
;; FONT-LOCK-KEYWORDS                                        |
;; 11/10/2006                                                |
;;-----------------------------------------------------------+
(defconst apollo-config-font-lock-keywords
  (list
   ;; The regexp-opt function optimizes the regular expression 
   ;; for a list of characters 

   ;;-----------------------------+
   ;; VARIABLE NAMES              |
   ;;-----------------------------+
   ;; In Apollo these are refered to as keys
   (cons (concat "\\<\\(" (regexp-opt '
			   ("Memory" 
			    "DataAdapterInstall"
			    "PublicSeqDbURL"
			    "PFetchServer"
			    "BrowserProgram"
			    "AskConfirmOverwrite"
			    "SiteShowLimit"
			    "FastDrawLimit"
			    "TextAvoidLimit"
			    "FrameOrientation"
			    "MainWindowWidth"
			    "Karyotypes"
			    "OutputTransactionXML"
			    "OutputChadoTransaction"
			    "ChadoTransactionTemplate"
			    "ChadoJdbcAdapterConfigFile"
			    "CommandLineXmlFileFormat"
			    "SpeciesToStyle"
			    "DO-ONE-LEVEL-ANNOTS"
			    "AdapterHistoryFile"
			    "AutosaveFile"
			    "JavaPath"
			    "EnableHttpServer"
			    "AutosaveInterval") 
			   t) "\\)\\>")'font-lock-variable-name-face)


   )
"Minimal highlighting expression for apollo config mode")

;;-----------------------------------------------------------+
;; COMMENTS IN APOLLO STYLES FILES                           |
;;-----------------------------------------------------------+
(defvar apollo-config-mode-syntax-table
  (let (( apollo-config-mode-syntax-table (make-syntax-table)))

    ;; < is the open comment
    ;; > is the close comment
    ;; The following thinks that / is a comment
    ;;(modify-syntax-entry ?/ "<" apollo-config-mode-syntax-table)
    (modify-syntax-entry ?/ ". 12" apollo-config-mode-syntax-table)
    (modify-syntax-entry ?\n ">" apollo-config-mode-syntax-table)
    apollo-config-mode-syntax-table)
  "Syntax table for apollo-config-mode")

;;-----------------------------------------------------------+
;; APOLLO TIERS MODE FUNCTION                                |
;;-----------------------------------------------------------+
;; This is the function that will be called by Emacs when 
;; the mode is started
(defun apollo-config-mode ()
  "Major mode for editing apollo config files"
  (interactive)
  (kill-all-local-variables)
  (set-syntax-table apollo-config-mode-syntax-table)
  (set (make-local-variable 'font-lock-defaults) 
       '(apollo-config-font-lock-keywords))
  (setq mode-name "Apollo Style")
  (run-hooks 'apollo-config-mode-hook)
)

(provide 'apollo-config-mode)