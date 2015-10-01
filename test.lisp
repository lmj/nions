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

(defpackage #:nions-test
  (:export #:run)
  (:use #:cl #:nions #:nions-test.1am))

(in-package :nions-test)

;;;; define algebras

(define-algebra
  :name my-complex
  :scalar-type fixnum
  :coords (c-real c-imag)
  :algebra-ops (c* c+ c- c-neg c-conj c-purepart))

(define-algebra/vector
  :power 1
  :scalar-ops (* + - -)
  :algebra-ops (vc* vc+ vc- vc-neg vc-conj vc-purepart))

;;; quaternions
(define-algebra
  :name quat
  :scalar-type integer
  :scalar-ops (* + - -)
  :coords (quat-r quat-i quat-j quat-k)
  :algebra-ops (quat* quat+ quat- quat-neg quat-conj quat-purepart))

(define-algebra/vector
  :power 2
  :scalar-type fixnum
  :scalar-ops (* + - -)
  :algebra-ops (vquat* vquat+ vquat- vquat-neg vquat-conj vquat-purepart))

;;; octonions
(define-algebra
  :name oct
  :scalar-type integer
  :scalar-ops (* + - -)
  :coords (oct0 oct1 oct2 oct3 oct4 oct5 oct6 oct7)
  :algebra-ops (oct* oct+ oct- oct-neg oct-conj oct-purepart))

;;;; manually defined operations

(defun manual-quat* (aw ax ay az bw bx by bz)
  (vector (+ (* aw bw) (- (* ax bx)) (- (* ay by)) (- (* az bz)))
          (+ (* aw bx)    (* ax bw)     (* ay bz)  (- (* az by)))
          (+ (* aw by) (- (* ax bz))    (* ay bw)     (* az bx))
          (+ (* aw bz)    (* ax by)  (- (* ay bx))    (* az bw))))

(defun manual-oct-conj (x)
  (oct (oct0 x)
       (- (oct1 x))
       (- (oct2 x))
       (- (oct3 x))
       (- (oct4 x))
       (- (oct5 x))
       (- (oct6 x))
       (- (oct7 x))))

(defun manual-oct-dot (x y)
  (+ (* (oct0 x) (oct0 y))
     (* (oct1 x) (oct1 y))
     (* (oct2 x) (oct2 y))
     (* (oct3 x) (oct3 y))
     (* (oct4 x) (oct4 y))
     (* (oct5 x) (oct5 y))
     (* (oct6 x) (oct6 y))
     (* (oct7 x) (oct7 y))))

;;;; util

(defun random-integers (count end)
  (loop repeat count
        collect (- (random (* 2 end)) end)))

(defun random-algebra (constructor dim)
  (apply constructor (random-integers dim 100)))

(defun vector/type (type &rest elems)
  (make-array (length elems) :element-type type :initial-contents elems))

(defmacro rtest (name &body body)
  `(test ,name
     (loop repeat 100
           do (progn ,@body))))

;;;; tests

(rtest complex-test
  (let* ((x (random-algebra #'my-complex 2))
         (y (random-algebra #'my-complex 2))
         (expected (* (complex (c-real x) (c-imag x))
                      (complex (c-real y) (c-imag y)))))
    (is (equalp (c-conj x)
                (my-complex (c-real x) (- (c-imag x)))))
    (is (eql expected
             (let ((z (c* x y)))
               (complex (c-real z) (c-imag z)))))
    (is (eql expected
             (let ((z (vc* (vector (c-real x) (c-imag x))
                           (vector (c-real y) (c-imag y)))))
               (complex (aref z 0) (aref z 1)))))))

(rtest quat-test
  (let* ((x (random-algebra #'quat 4))
         (y (random-algebra #'quat 4))
         (expected (manual-quat* (quat-r x) (quat-i x) (quat-j x) (quat-k x)
                                 (quat-r y) (quat-i y) (quat-j y) (quat-k y))))
    (is (equalp (quat-conj x)
                (quat (quat-r x)
                      (- (quat-i x))
                      (- (quat-j x))
                      (- (quat-k x)))))
    (is (equalp expected
                (let ((z (quat* x y)))
                  (vector (quat-r z)
                          (quat-i z)
                          (quat-j z)
                          (quat-k z)))))
    (is (equalp expected
                (vquat* (vector/type 'fixnum
                                     (quat-r x)
                                     (quat-i x)
                                     (quat-j x)
                                     (quat-k x))
                        (vector/type 'fixnum
                                     (quat-r y)
                                     (quat-i y)
                                     (quat-j y)
                                     (quat-k y)))))))

(rtest oct-test
  (let ((x (random-algebra #'oct 8))
        (y (random-algebra #'oct 8)))
    (is (equalp (manual-oct-conj x)
                (oct-conj x)))
    (is (zerop (oct0 (oct- (oct* y (oct-conj x))
                           (oct* x (oct-conj y))))))
    (let ((z (oct+ (oct* y (oct-conj x))
                   (oct* x (oct-conj y)))))
      (is (= (manual-oct-dot x y)
             (/ (oct0 z) 2)))
      (is (zerop (oct1 z)))
      (is (zerop (oct2 z)))
      (is (zerop (oct3 z)))
      (is (zerop (oct4 z)))
      (is (zerop (oct5 z)))
      (is (zerop (oct6 z)))
      (is (zerop (oct7 z))))))

(test error-test
  (signals error
    (macroexpand-1 '(define-algebra
                     :name scalar
                     :coords (value)
                     :scalar-type fixnum
                     :algebra-ops (s* s+ s- s-neg s-conj s-purepart))))
  (signals error
    (macroexpand-1 '(define-algebra/vector
                     :power 0
                     :scalar-type fixnum
                     :algebra-ops (s* s+ s- s-neg s-conj s-purepart))))
  (signals error
    (macroexpand-1 '(define-algebra
                     :name tri
                     :coords (p q r)
                     :scalar-type fixnum
                     :algebra-ops (s* s+ s- s-neg s-conj s-purepart)))))
