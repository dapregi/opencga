package org.opencb.opencga.analysis.variant.manager;

import org.apache.commons.lang3.StringUtils;
import org.opencb.commons.datastore.core.ObjectMap;
import org.opencb.opencga.analysis.ConfigurationUtils;
import org.opencb.opencga.core.exception.ToolExecutorException;
import org.opencb.opencga.catalog.exceptions.CatalogException;
import org.opencb.opencga.catalog.managers.CatalogManager;
import org.opencb.opencga.core.config.Configuration;
import org.opencb.opencga.storage.core.StorageEngineFactory;
import org.opencb.opencga.storage.core.config.StorageConfiguration;
import org.opencb.opencga.storage.core.exceptions.StorageEngineException;
import org.opencb.opencga.storage.core.variant.VariantStorageEngine;

import java.io.IOException;

/**
 * Helper interface to be used by opencga local analysis executors.
 */
public interface VariantStorageAnalysisExecutor {

    ObjectMap getExecutorParams();


    default VariantStorageManager getVariantStorageManager() throws ToolExecutorException {
        String opencgaHome = getExecutorParams().getString("opencgaHome");
        try {
            Configuration configuration = ConfigurationUtils.loadConfiguration(opencgaHome);
            StorageConfiguration storageConfiguration = ConfigurationUtils.loadStorageConfiguration(opencgaHome);

            CatalogManager catalogManager = new CatalogManager(configuration);
            StorageEngineFactory engineFactory = StorageEngineFactory.get(storageConfiguration);
            return new VariantStorageManager(catalogManager, engineFactory);
        } catch (CatalogException | IOException e) {
            throw new ToolExecutorException(e);
        }
    }

    default VariantStorageEngine getVariantStorageEngine() throws ToolExecutorException {
        ObjectMap executorParams = getExecutorParams();
        String storageEngine = executorParams.getString("storageEngineId");
        String dbName = executorParams.getString("dbName");
        if (StringUtils.isEmpty(storageEngine) || StringUtils.isEmpty(dbName)) {
            throw new ToolExecutorException("Missing arguments!");
        } else {
            try {
                return StorageEngineFactory.get().getVariantStorageEngine(storageEngine, dbName);
            } catch (StorageEngineException e) {
                throw new ToolExecutorException(e);
            }

        }
    }

}
