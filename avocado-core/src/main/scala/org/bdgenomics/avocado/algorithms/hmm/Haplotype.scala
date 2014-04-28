package org.bdgenomics.avocado.algorithms.hmm

import scala.math.Ordering
import org.bdgenomics.adam.rich.RichADAMRecord
import scala.math._

/**
 * Haplotype generated from HMM alignment.
 *
 * @param sequence String representing haplotype alignment.
 */
class Haplotype(val sequence: String, region: Seq[RichADAMRecord], hmm: HMMAligner = new HMMAligner, val reference: String = "") {

  lazy val referenceAlignment = hmm.alignSequences(reference, sequence, null)
  lazy val hasVariants = referenceAlignment.hasVariants

  lazy val perReadLikelihoods: Seq[Double] = region.map(read => {
    try {
      val alignment = HMMAligner.align(sequence, read.getSequence.toString, null)
      alignment.likelihood + alignment.prior
    } catch {
      case _: Throwable => {
        0.0
      }
    }
  })

  lazy val readsLikelihood = perReadLikelihoods.sum

  override def toString(): String = {
    sequence
  }

}

/**
 * Haplotypes are ordered by increasing reads likelihood, assuming they
 * come from the same group of reads.
 */
object HaplotypeOrdering extends Ordering[Haplotype] {

  /**
   * Compares two haplotypes. Returns (-1, 0, 1) if h1 has (lower, same, higher) read
   * likelihood than h2. 0 is only returned if they have the same sequence AND likelihood
   *
   * @param h1 First haplotype to compare.
   * @param h2 Second haplotype to compare.
   * @return Ordering info for haplotypes.
   */
  def compare(h1: Haplotype, h2: Haplotype): Int = {
    if (h1.sequence == h2.sequence) {
      h1.readsLikelihood.compare(h2.readsLikelihood)
    } else if (h1.readsLikelihood < h2.readsLikelihood) {
      -1
    } else {
      1
    }
  }
}

