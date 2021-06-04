(require 'package)
(defvar my-packages '(ox-rst org-ref gnuplot))
(add-to-list 'package-archives
	     '("melpa" . "https://melpa.org/packages/"))

(package-initialize)

(if (version< emacs-version "26.3")
  (progn
    (message "Emacs < 26.3, setting gnutls-algorithm-priority TLS1.3 ")
    (setq gnutls-algorithm-priority "NORMAL:-VERS-TLS1.3"))
  (progn
    (message "Emacs >= 26.3, not setting gnutls-algorithm-priority TLS1.3")
    ;(setq gnutls-algorithm-priority "NORMAL:-VERS-TLS1.3")
    )
)

(when (not package-archive-contents)
  (package-refresh-contents))

(dolist (p my-packages)
  (when (not (package-installed-p p))
    (package-install p)))

;; forceful package-refresh
(dolist (p my-packages)
  (when (not (package-installed-p p))
    (package-refresh-contents)
    (package-install p)))

;; (source: http://emacs.stackexchange.com/a/22174/93)
(defun package-update (package &optional version)
  "Update a package from the package archives.
If VERSION is nil, update the package to the most recent version
available.  Otherwise, VERSION should be a version string, or a
list of the type returned by `version-to-list'. The package will
be updated only if the currently installed version is less than
the version specified, even if a newer version is available."
  (unless package--initialized
    (package-initialize t))
  (unless package-archive-contents
    (package-refresh-contents))
  (let* ((current (cadr (assoc package package-alist)))
         (current-version (if current (package-desc-version current) '(-1)))
         (pkg-desc (cadr (assoc package package-archive-contents)))
         (pkg-version (and pkg-desc (package-desc-version pkg-desc)))
         (target-version (or (and (stringp version) (version-to-list version))
                             version
                             pkg-version)))
    (when (version-list-< current-version target-version)
      (if (null pkg-desc)
        (error "Package `%s' not found in package archives" package))
      (if (version-list-< pkg-version target-version)
        (error "A suitable version of `%s' is not available" package))
      (package-install pkg-desc)
      (if current
          (package-delete current)))
    nil))

(package-update 'org) ; '(9 4 0))

(provide 'emacs-install-pkgs)
