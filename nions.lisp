;;; Copyright (c) 2015 James M. Lawrence
;;;
;;; Permission is hereby granted, free of charge, to any person
;;; obtaining a copy of this software and associated documentation
;;; files (the "Software"), to deal in the Software without
;;; restriction, including without limitation the rights to use, copy,
;;; modify, merge, publish, distribute, sublicense, and/or sell copies
;;; of the Software, and to permit persons to whom the Software is
;;; furnished to do so, subject to the following conditions:
;;;
;;; The above copyright notice and this permission notice shall be
;;; included in all copies or substantial portions of the Software.
;;;
;;; THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
;;; EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
;;; MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
;;; NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
;;; HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
;;; WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
;;; OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
;;; DEALINGS IN THE SOFTWARE.

(defpackage #:nions
  (:export #:define-algebra
           #:define-algebra/vector)
  (:use #:cl))

(in-package :nions)

;;;; scalar operations

(defstruct (scalar-ops (:constructor scalar-ops (* + - neg))
                       (:conc-name scalar-))
  * + - neg)

;;; expansion-time-only environment; contains a scalar-ops instance
(defvar *scalar-ops*)

;;;; expanders -- each with 'e' prefix

(defun %e-left/right (form symbol pos)
  (if (and (consp form) (eq (first form) 'pair))
      (funcall pos form)
      `(,symbol ,form)))

(defun e-left (form)
  "Expand left part of algebra element."
  (%e-left/right form 'L #'second))

(defun e-right (form)
  "Expand right part of algebra element."
  (%e-left/right form 'R #'third))

(defun e-neg (a n)
  "Expand negation operation."
  (if (zerop n)
      `(,(scalar-neg *scalar-ops*) ,a)
      (let ((k (- n 1)))
        `(pair ,(e-neg (e-left  a) k)
               ,(e-neg (e-right a) k)))))

(defun e-conj (a n)
  "Expand conjugate operation."
  (if (zerop n)
      a
      (let ((k (- n 1)))
       `(pair ,(e-conj (e-left  a) k)
              ,(e-neg  (e-right a) k)))))

(defun %e-+/- (a b n scalar-op)
  (if (zerop n)
      `(,scalar-op ,a ,b)
      (let ((k (- n 1)))
        `(pair ,(%e-+/- (e-left  a) (e-left  b) k scalar-op)
               ,(%e-+/- (e-right a) (e-right b) k scalar-op)))))

(defun e+ (a b n)
  "Expand + operation."
  (%e-+/- a b n (scalar-+ *scalar-ops*)))

(defun e- (a b n)
  "Expand - operation."
  (%e-+/- a b n (scalar-- *scalar-ops*)))

(defun e* (a b n)
  "Expand * operation."
  (if (zerop n)
      `(,(scalar-* *scalar-ops*) ,a ,b)
      (let ((k (- n 1)))
        `(pair ,(e- (e*         (e-left  a)    (e-left  b) k)
                    (e* (e-conj (e-right b) k) (e-right a) k)
                    k)
               ,(e+ (e* (e-right b)         (e-left a)    k)
                    (e* (e-right a) (e-conj (e-left b) k) k)
                    k)))))

(defun e-identity (a n)
  "Expand identity operation."
  (if (zerop n)
      a
      (let ((k (- n 1)))
        `(pair ,(e-identity (e-left a) k)
               ,(e-identity (e-right a) k)))))

(defun e-purepart (a n)
  "Expand purepart operation."
  (if (zerop n)
      0
      (let ((k (- n 1)))
        `(pair ,(e-purepart (e-left a) k)
               ,(e-identity (e-right a) k)))))

;;;; convert from LR coordinates to regular coordinates

(defun binary-digits->integer (digits)
  "(0 0 1 1) --> 12"
  (loop with result = 0
        for digit in digits
        for n from 0
        do (ecase digit
             (0)
             (1 (incf result (expt 2 n))))
        finally (return result)))

(defun pairs->list (code)
  "Produce a flat list from recursive pairs. 
  (pair (pair a b) (pair c d)) --> (a b c d)"
  (if (and (consp code) (eq (first code) 'pair))
      (nconc (pairs->list (second code))
             (pairs->list (third code)))
      (list code)))

(defun extract-lrs (code)
  "(L (R (L (L X)))) --> (L R L L), X"
  (loop for pos = code then (second pos)
        while (consp pos)
        collect (first pos) into lrs
        finally (return (values lrs pos))))

(defun lrs->integer (lrs)
  "L and R are binary digits 0 and 1, respectively. 
   (L R R L) --> 6"
  (binary-digits->integer (mapcar (lambda (x)
                                    (ecase x
                                      (L 0)
                                      (R 1)))
                                  lrs)))

(defun convert-expr (emit-coord expr)
  "Convert all LR forms to coordinate forms.
   (* (L (R (R (L X)))) (L (L (L (L Y))))) --> (* (AREF X 6) (AREF Y 0))"
  (if (consp expr)
      (destructuring-bind (first . rest) expr
        (if (member first '(L R))
            (multiple-value-bind (lrs sym) (extract-lrs expr)
              (funcall emit-coord sym (lrs->integer lrs)))
            (cons (convert-expr emit-coord first)
                  (convert-expr emit-coord rest))))
      expr))

(defun convert-alg (alg emit-alg emit-coord)
  "Convert an entire algebra element to coordinate form."
  (funcall emit-alg (mapcar (lambda (expr) (convert-expr emit-coord expr))
                            (pairs->list alg))))

;;;; generate defuns

(defmacro free-symbol (symbol)
  "For cleaner and more compact expansions via avoiding gensyms."
  `(make-symbol ,(symbol-name symbol)))

(defun def-op (name expander arity power emit-alg emit-coord type)
  "Generate a `defun' for the given expander."
  (let ((args (ecase arity
                (1 `(,(free-symbol a)))
                (2 `(,(free-symbol a) ,(free-symbol b))))))
    `(defun ,name ,args
       (declare (type ,type ,@args))
       ,(convert-alg (apply expander (append args (list power)))
                     emit-alg emit-coord))))

(defun def-ops (power scalar-ops algebra-ops emit-alg emit-coord type)
  "Generate all `defun's."
  (destructuring-bind (alg* alg+ alg- alg-neg alg-conj alg-purepart) algebra-ops
    (let ((*scalar-ops* (apply #'scalar-ops scalar-ops)))
      `(progn
         ,(def-op alg*         'e*         2 power emit-alg emit-coord type)
         ,(def-op alg+         'e+         2 power emit-alg emit-coord type)
         ,(def-op alg-         'e-         2 power emit-alg emit-coord type)
         ,(def-op alg-neg      'e-neg      1 power emit-alg emit-coord type)
         ,(def-op alg-conj     'e-conj     1 power emit-alg emit-coord type)
         ,(def-op alg-purepart 'e-purepart 1 power emit-alg emit-coord type)))))

;;;; front end

(defun power-of-two-p (n)
  (and (plusp n)
       (loop for k below (- (integer-length n) 1)
             never (logbitp k n))))

(defun check-dim (dim definer)
  (unless (> dim 1)
    (error "In ~a: dimension must be greater than 1. " definer))
  (unless (power-of-two-p dim)
    (error "In ~a: dimension is not a power of two." definer)))

(defun def-print-object (name coords)
  `(defmethod print-object ((alg ,name) stream)
     (let ((*package* (find-package 'keyword)))
       (format stream "(~s ~{~s~^ ~})" ',name
               (list ,@(loop for coord in coords
                             collect `(,coord alg)))))))

(defmacro define-algebra (&key name coords algebra-ops
                               (scalar-ops '(* + - -)) 
                               (scalar-type t))
  "Define a Cayley-Dickson algebra.

  `name' -- Name of the algebra. This will be the constructor.

  `coords' -- A list of accessors to be defined by this macro, giving
  coordinate values (the coefficients of the basis elements). The
  number of coordinates must be a power of two.

  `algebra-ops' -- A list of algebra operations to be defined by this
  macro, respectively: multiply, add, subtract, negate, conjugate,
  purepart.

  `scalar-ops' -- A list of operations to be used for scalars:
  multiply, add, subtract, negate. For built-in scalars the default
  value sufficies: (* + - -).

  `scalar-type' -- Optional type for scalars, such as fixnum."
  (flet ((emit-alg (forms)
           `(,name ,@(loop for form in forms
                           collect `(the ,scalar-type ,form))))
         (emit-coord (alg index)
           `(,(elt coords index) ,alg)))
    (let ((dim (length coords)))
      (check-dim dim 'define-algebra)
      (let ((power (round (log dim 2))))
        `(progn
           (defstruct (,name (:constructor ,name ,coords)
                             (:conc-name nil))
             ,@(loop for coord in coords
                     collect `(,coord nil :type ,scalar-type)))
           ,(def-ops power scalar-ops algebra-ops #'emit-alg #'emit-coord name)
           ,(def-print-object name coords))))))

(defmacro define-algebra/vector (&key power algebra-ops
                                      (scalar-ops '(* + - -)) 
                                      (scalar-type t))
  "Define a Cayley-Dickson algebra with vector representation.

  `power' -- The dimension of the algebra will be (expt 2 power).

  `algebra-ops' -- A list of algebra operations to be defined by this
  macro, respectively: multiply, add, subtract, negate, conjugate,
  purepart.

  `scalar-ops' -- A list of operations to be used for scalars:
  multiply, add, subtract, negate. For built-in scalars the default
  value sufficies: (* + - -).

  `scalar-type' -- Optional type for scalars, such as fixnum."
  (let ((dim (expt 2 power)))
    (check-dim dim 'define-algebra/vector)
    (flet ((emit-alg (forms)
             (let ((result (free-symbol result)))
               `(let ((,result (make-array ,dim :element-type ',scalar-type)))
                  ,@(loop for form in forms
                          for index from 0
                          collect `(setf (aref ,result ,index)
                                         (the ,scalar-type ,form)))
                  ,result)))
           (emit-coord (alg index)
             `(aref ,alg ,index)))
      (def-ops power scalar-ops algebra-ops #'emit-alg #'emit-coord
               `(simple-array ,scalar-type (,dim)))))) 
