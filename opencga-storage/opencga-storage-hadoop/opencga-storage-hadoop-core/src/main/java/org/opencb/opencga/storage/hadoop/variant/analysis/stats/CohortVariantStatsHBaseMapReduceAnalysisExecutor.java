package org.opencb.opencga.storage.hadoop.variant.analysis.stats;

import org.opencb.commons.datastore.core.ObjectMap;
import org.opencb.commons.datastore.core.Query;
import org.opencb.commons.utils.StringUtils;
import org.opencb.opencga.core.tools.variant.CohortVariantStatsAnalysisExecutor;
import org.opencb.opencga.core.annotations.ToolExecutor;
import org.opencb.opencga.core.common.TimeUtils;
import org.opencb.opencga.core.exception.ToolException;
import org.opencb.opencga.core.exception.ToolExecutorException;
import org.opencb.opencga.storage.core.exceptions.StorageEngineException;
import org.opencb.opencga.storage.core.metadata.models.TaskMetadata;
import org.opencb.opencga.storage.core.variant.adaptors.VariantQueryException;
import org.opencb.opencga.storage.hadoop.variant.HadoopVariantStorageEngine;
import org.opencb.opencga.storage.hadoop.variant.adaptors.VariantHadoopDBAdaptor;
import org.opencb.opencga.storage.hadoop.variant.analysis.HadoopVariantAnalysisExecutor;
import org.opencb.opencga.storage.hadoop.variant.stats.CohortVariantStatsDriver;

@ToolExecutor(id = "hbase-mapreduce", tool = "cohort-variant-stats",
        framework = ToolExecutor.Framework.MAP_REDUCE,
        source = ToolExecutor.Source.HBASE)
public class CohortVariantStatsHBaseMapReduceAnalysisExecutor
        extends CohortVariantStatsAnalysisExecutor implements HadoopVariantAnalysisExecutor {

    @Override
    public void run() throws ToolException {

        HadoopVariantStorageEngine engine = getHadoopVariantStorageEngine();

        VariantHadoopDBAdaptor dbAdaptor;
        int studyId;
        Integer cohortId;

        String temporaryCohortName = "TEMP_" + StringUtils.randomString(5) + "_" + TimeUtils.getTimeMillis();
        try {
            dbAdaptor = engine.getDBAdaptor();

            studyId = dbAdaptor.getMetadataManager().getStudyId(getStudy());

            cohortId = dbAdaptor.getMetadataManager().registerCohort(getStudy(), temporaryCohortName, getSampleNames());
            dbAdaptor.getMetadataManager().updateCohortMetadata(studyId, cohortId, cohortMetadata -> {
                cohortMetadata.setStatus("TEMPORARY", TaskMetadata.Status.RUNNING);
                return cohortMetadata;
            });
        } catch (StorageEngineException e) {
            throw new ToolExecutorException(e);
        }

        try {
//            Query variantsQuery = getVariantsQuery();
            Query variantsQuery = new Query();

            ObjectMap params = new ObjectMap(variantsQuery)
                    .append(CohortVariantStatsDriver.COHORT, temporaryCohortName)
                    .append(CohortVariantStatsDriver.OUTPUT, getOutputFile().toAbsolutePath().toUri());
            engine.getMRExecutor().run(CohortVariantStatsDriver.class, CohortVariantStatsDriver.buildArgs(
                    dbAdaptor.getArchiveTableName(studyId),
                    dbAdaptor.getVariantTable(),
                    studyId,
                    null,
                    params
            ), engine.getOptions(), "Calculate cohort variant stats");

        } catch (VariantQueryException | StorageEngineException e) {
            throw new ToolExecutorException(e);
        } finally {
            dbAdaptor.getMetadataManager().removeCohort(studyId, cohortId);
        }
    }

}
