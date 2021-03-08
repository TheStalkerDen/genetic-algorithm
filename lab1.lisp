(uiop:define-package :lab1-genetic-algo
  (:nicknames :lab1))

(in-package :lab1-genetic-algo)

(declaim (optimize (safety 3)))

(defparameter *arg-max* 10)

(defvar *max-limit*)
(defvar *min-limit*)
(defvar *current-fun*)
(defvar *tournament-chroms-count*)
(defvar *elit-elements-count*)
(defvar *mutation-chance*)

;;;all about chromosomes
(defclass chromosome ()
  ((genes
    :initarg :genes
    :accessor genes)
   (z-val
    :initarg :z-val
    :initform 0
    :accessor z-val
    :type number
    :documentation "z = f(genes)")))

(defmethod print-object ((obj chromosome) stream)
  (with-slots (genes z-val) obj
    (print-unreadable-object
        (obj stream)
      (format stream "Chromosome z=~A ~A" z-val genes))))

(defun generate-genes ()
  (loop :repeat *arg-max*
        :collect (+ (random (abs (- *min-limit*
                                    *max-limit*)))
                    *min-limit*)
          :into res
        :finally (return (coerce res 'vector))))

(defun make-chromosome (&key genes)
  (make-instance
   'chromosome
   :genes (if genes
              genes
              (generate-genes))))

;;;convenient functions for testing
(defun test-sphere-function (population-size
                             iterations
                             &key (elitism-percentage 5) (mutation-chance 5))
  (genetic-algo-start #'sphere-function population-size
                      -100.0
                      100.0
                      iterations
                      :elitism-percentage elitism-percentage
                      :mutation-chance mutation-chance))

(defun test-ackley-function (population-size
                             iterations
                             &key (elitism-percentage 5) (mutation-chance 5))
 (genetic-algo-start #'ackley-function population-size
                      -32.768
                      32.768
                      iterations
                      :elitism-percentage elitism-percentage
                      :mutation-chance mutation-chance))

(defun test-griewank-function (population-size
                               iterations
                               &key (elitism-percentage 5) (mutation-chance 5))
 (genetic-algo-start #'griewank-function population-size
                      -600.0
                      600.0
                      iterations
                      :elitism-percentage elitism-percentage
                      :mutation-chance mutation-chance))

(defun test-rastrigin-function (population-size
                                iterations
                                &key (elitism-percentage 5) (mutation-chance 5))
 (genetic-algo-start #'rastrigin-function population-size
                      -5.12
                      5.12
                      iterations
                      :elitism-percentage elitism-percentage
                      :mutation-chance mutation-chance))

(defun test-rosenbrock-function (population-size
                                 iterations
                                 &key (elitism-percentage 5) (mutation-chance 5))
 (genetic-algo-start #'rosenbrock-function population-size
                      -5.0
                      10.0
                      iterations
                      :elitism-percentage elitism-percentage
                      :mutation-chance mutation-chance))

;;;main functions of genetic algorithm
(defun genetic-algo-start (fun population-size min-limit max-limit iterations
                           &key (elitism-percentage 5) (mutation-chance 5))
  (let* ((*max-limit* max-limit)
         (*min-limit* min-limit)
         (*current-fun* fun)
         (*elit-elements-count* (ceiling (* population-size
                                            (/ elitism-percentage
                                               100))))
         (*mutation-chance* mutation-chance)
         (*tournament-chroms-count* (- population-size
                                       *elit-elements-count*))
         (start-population
           (create-start-population population-size))
         (final-population
           (loop :repeat iterations
                 :for population = (population-estimation start-population)
                   :then (population-estimation population)
                 :do (setf population (middle-processes population))
                 :finally (return (population-estimation population)))))
    (first (sort-population final-population))))

(defun create-start-population (n)
  (loop :repeat n
        :collect (make-chromosome)))

(defun sort-population (population)
  (sort population (lambda (z-val1 z-val2)
                     (< z-val1 z-val2))
        :key (lambda (chrom)
               (z-val chrom))))

(defun population-estimation (population)
  (dolist (chrom population)
    (setf (z-val chrom) (funcall *current-fun* (genes chrom))))
  population)

(defun middle-processes (population)
  (multiple-value-bind (elits tournament-chroms)
      (selection population)
    `(,@elits ,@(mutation (flat-crossover tournament-chroms)))))

(defun selection (population)
  (let ((sorted-chromes (sort-population population)))
    (values
     (subseq sorted-chromes 0 *elit-elements-count*)
     (tournament-hunger-game (nthcdr *elit-elements-count* sorted-chromes)))))

;;;this is a joke-name for tournament selection :)
(defun tournament-hunger-game (tournament-participants)
  (loop :repeat *tournament-chroms-count*
        :with tourn-vector = (apply #'vector tournament-participants)
        :collect (let ((chrom1 (aref tourn-vector
                                     (random *tournament-chroms-count*)))
                       (chrom2 (aref tourn-vector
                                     (random *tournament-chroms-count*))))
                   (if (< (z-val chrom1) (z-val chrom2))
                       chrom1
                       chrom2))))
  
(defun flat-crossover (tournament-chroms)
  (loop :repeat *tournament-chroms-count*
        :with tourn-vector = (apply #'vector tournament-chroms)
        :collect (let ((chrom1 (aref tourn-vector
                                     (random *tournament-chroms-count*)))
                       (chrom2 (aref tourn-vector
                                     (random *tournament-chroms-count*))))
                   (make-chromosome
                    :genes (loop
                             :for c1 :across (genes chrom1)
                             :and c2 :across (genes chrom2)
                             :collect (if (= (- c1 c2) 0)
                                          c1 
                                          (+ (random (abs (- c1 c2)))
                                             (min c1 c2)))
                               :into lst
                             :finally (return (coerce lst 'vector)))))))

(defun mutation (crossover-results)
  (dolist (chrom crossover-results)
    (when (< (random 100) *mutation-chance*)
      (setf (aref (genes chrom) (random *arg-max*))
            (+ (random (abs (- *max-limit*
                               *min-limit*)))
               *min-limit*))))
  crossover-results)


(defun pow (val degree)
  (if (= degree 0)
      1
      (* val (pow val (1- degree)))))

;;;benchmark functions

(defun sphere-function (val-vector)
  (loop :for x :across val-vector
        :summing (pow x 2)))

(defun ackley-function (val-vector)
  (let ((part1
          (/ (loop :for x :across val-vector
                   :summing (pow x 2))
             (length val-vector)))
        (part2
          (/ (loop :for x :across val-vector
                   :summing (cos (* 2 pi x)))
             (length val-vector))))
    (+ (- (* -20
             (exp (* -0.2
                     (sqrt part1))))
          (exp part2))
       20
       (exp 1))))

(defun griewank-function (val-vector)
  (let ((part1
          (/ (loop :for x :across val-vector
                   :summing (pow x 2))
             4000))
        (part2
          (loop :for x :across val-vector
                :and i = 1 :then (1+ i)
                :for res = 1
                  :then (* res
                           (cos (/ x (sqrt i))))
                :finally (return res))))
    (+ (- part1 part2) 1)))

(defun rastrigin-function (val-vector)
  (let ((part1
          (loop :for x :across val-vector
                :summing (- (pow x 2)
                            (* 10
                               (cos (* 2 pi x)))))))
    (+ (* 10
          (length val-vector))
       part1)))

(defun rosenbrock-function (val-vector)
  (loop :for i :from 1 :below (length val-vector)
        :summing (+ (* 100
                       (pow (- (aref val-vector i)
                               (pow
                                (aref val-vector (1- i))
                                2))
                            2))
                    (pow (1- (aref val-vector (1- i)))
                         2))))
