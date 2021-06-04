(defun my/return-t (orig-fun &rest args)
t)
(defun my/disable-yornp (orig-fun &rest args)
  (advice-add 'yes-or-no-p :around #'my/return-t)
  (advice-add 'y-or-n-p :around #'my/return-t)
  (let ((res (apply orig-fun args)))
    (advice-remove 'yes-or-no-p #'my/return-t)
    (advice-remove 'y-or-n-p #'my/return-t)
    res))
(advice-add 'org-babel-tangle :around #'my/disable-yornp)

;(setq out-dir "/tmp/orgie2")


(defun org-export-output-file-name-modified
    (orig-fun extension &optional subtreep pub-dir)
  (unless pub-dir
    (setq pub-dir out-dir)
    (unless (file-directory-p pub-dir)
      (make-directory pub-dir)))
  (apply orig-fun extension subtreep pub-dir nil))

(advice-add 'org-export-output-file-name
            :around #'org-export-output-file-name-modified)


(package-initialize)

;(add-to-list 'load-path ".")
;(require 'emacs-install-pkgs)

(if (eq system-type 'windows-nt)
  (setq org-plantuml-jar-path "C:/ProgramData/chocolatey/lib/plantuml/tools/plantuml.jar")
  (setq org-plantuml-jar-path "/usr/share/plantuml/plantuml.jar")
)

(setq plantuml-default-exec-mode 'jar)
(setq my-languages '((shell . t)
                     (python . t)
                     (gnuplot . t)
                     (plantuml . t)
                     (C . t)))
(org-babel-do-load-languages 'org-babel-load-languages my-languages)
(setq auto-save-default nil make-backup-files nil)

;; Disable a lot of export parameters,https://orgmode.org/manual/Export-Settings.html
(setq org-export-with-broken-links t)
(setq org-export-with-author nil)
(setq org-export-with-date nil)
(setq org-export-with-email nil)
(setq org-export-with-section-numbers nil)
(setq org-export-with-toc nil)

(setq org-rst-link-use-ref-role nil)
(setq org-rst-code-block 'code-block)
;;(setq org-ref-default-bibliography '("~/bibs/bib.bib"))

(require 'org)
(require 'org-ref)
(require 'ox-rst)

(defun my/org-ref-rst-export (orig-func keyword desc format)
  "Add rst export for ORIG-FUNC, use KEYWORD, DESC, FORMAT."
  (if (eq format 'rst)
    (format ":ref:`%s`" keyword)
    (apply orig-func keyword desc format nil)))

(advice-add 'org-ref-ref-export :around #'my/org-ref-rst-export)

(defun my/org-ref-eqref-rst-export (orig-func keyword desc format)
  "Add rst export for ORIG-FUNC, use KEYWORD, DESC, FORMAT."
  (if (eq format 'rst)
    (format ":eq:`%s`" keyword)
    (apply orig-func keyword desc format nil)))

(advice-add 'org-ref-eqref-export :around #'my/org-ref-eqref-rst-export)

(defun my/org-cref-rst-export (orig-func keyword desc format)
  "Add rst export for ORIG-FUNC, use KEYWORD, DESC, FORMAT."
  (if (eq format 'rst)
    (format ":numref:`%s`" keyword)
    (apply orig-func keyword desc format nil)))

(advice-add 'org-ref-cref-export :around #'my/org-cref-rst-export)

(defun my/org-Cref-rst-export (orig-func keyword desc format)
  "Add rst export for ORIG-FUNC, use KEYWORD, DESC, FORMAT."
  (if (eq format 'rst)
    (format ":numref:`%s`" keyword)
    (apply orig-func keyword desc format nil)))

(advice-add 'org-ref-Cref-export :around #'my/org-Cref-rst-export)

(defun my/org-ref-format-cite (orig-func keyword desc format)
  "Add rst export for ORIG-FUNC, use KEYWORD, DESC, FORMAT."
  (cond ((eq format 'rst)
          (concat ":cite:`"
            (mapconcat
              (lambda (key)
                (format "%s" key))
              (org-ref-split-and-strip-string keyword) ",") "`"))
    (t (apply orig-func keyword desc format nil))))

(advice-add 'org-ref-format-cite :around #'my/org-ref-format-cite)

(defun my/org-rst-latex-environment (orig-func latex-environment _contents info)
  "Override default ORIG-FUNC org-rst to wrap LATEX-ENVIRONMENT in math directive, use INFO, not _CONTENTS."
  (let* ((attr (org-export-read-attribute :attr_rst latex-environment))
         (env_data (org-element-property :value latex-environment))
         (label (org-element-property :name latex-environment))
         (nowrap (plist-get attr :nowrap)))
    (cond
      ((plist-get info :with-latex)
        (concat
          ".. math::\n"
          (when nowrap (format "%s\n" (org-rst--indent-string ":nowrap:" org-rst-quote-margin)))
          (when label (format "%s %s\n" (org-rst--indent-string ":label:" org-rst-quote-margin) label))
          "\n"
          (format "%s"
            (org-rst--indent-string
              (org-remove-indentation env_data) org-rst-quote-margin))))
      (t env_data))
    ))

(advice-add 'org-rst-latex-environment :around #'my/org-rst-latex-environment)

(defun my/org-rst-src-block (orig-func src-block _contents info)
  "Transcode a SRC-BLOCK element from Org to reStructuredText.
CONTENTS holds the contents of the item.  INFO is a plist holding
contextual information."
  (when (org-string-nw-p (org-element-property :value src-block))
    (let* ((lang (org-element-property :language src-block))
		   (label (org-element-property :name src-block))
		   (value (org-remove-indentation
				   (org-element-property :value src-block)))
           (num-start (org-export-get-loc src-block info))
           (codeblockd (plist-get info :rst-code-block))
            (caption (org-export-get-caption src-block))
		   (attributes
			(org-export-read-attribute :attr_rst src-block))
		   (class (plist-get attributes :class)))
      (cond
       ;; Case 1. We only override this implementation
        ((eq codeblockd 'code-block)
          (let ((lst-lang
                  (or (cadr (assq (intern lang) org-rst-pygments-langs)) lang)))
            (concat
              (format ".. code-block:: %s\n" lst-lang)
              (when num-start (format "    :lineno-start: %s\n" (1+ num-start)))
              (when class (format "    :class: %s\n" class))
              (when label (format "    :name: %s\n" label))
              (when caption (format "    :caption: %s\n" caption))
              "\n"
              (org-rst--indent-string value org-rst-quote-margin))))
        (t (apply orig-func src-block _contents info nil))))))

(advice-add 'org-rst-src-block :around #'my/org-rst-src-block)


(defun jemacs-export-org-doc ()
  (dolist (file command-line-args-left)
    (princ (concat "Exporting " file " to " out-dir "/\n") #'external-debugging-output)
    (with-current-buffer (find-file-noselect file)
      (org-rst-export-to-rst))))
(advice-add 'jemacs-export-org-doc :around #'my/disable-yornp)
