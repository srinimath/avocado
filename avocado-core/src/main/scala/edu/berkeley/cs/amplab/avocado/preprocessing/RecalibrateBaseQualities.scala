/*
 * Copyright (c) 2014. Regents of the University of California
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

package edu.berkeley.cs.amplab.avocado.preprocessing

import edu.berkeley.cs.amplab.adam.avro.ADAMRecord
import edu.berkeley.cs.amplab.adam.models.SnpTable
import edu.berkeley.cs.amplab.adam.rdd.AdamContext._
import java.io.File
import org.apache.commons.configuration.SubnodeConfiguration
import org.apache.spark.rdd.RDD

object RecalibrateBaseQualities extends PreprocessingStage {
  
  val stageName = "recalibrateBaseQualities"

  def apply (rdd: RDD[ADAMRecord], config: SubnodeConfiguration): RDD[ADAMRecord] = {
    // check for snp table
    val snpTable = if (config.containsKey("snpTable")) {
      SnpTable(new File(config.getString("snpTable")))
    } else {
      SnpTable()
    }
    
    // run bqsr with snp table loaded
    rdd.adamBQSR(snpTable)
  }

}
