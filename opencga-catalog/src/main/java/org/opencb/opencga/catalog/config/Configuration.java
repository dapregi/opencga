/*
 * Copyright 2015-2016 OpenCB
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

package org.opencb.opencga.catalog.config;

import com.fasterxml.jackson.databind.ObjectMapper;
import com.fasterxml.jackson.dataformat.yaml.YAMLFactory;
import org.opencb.opencga.catalog.models.Project;
import org.opencb.opencga.catalog.models.acls.permissions.StudyAclEntry;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.util.List;

/**
 * Created by imedina on 16/03/16.
 */
public class Configuration {

    private String logLevel;
    private String logFile;

    private boolean openRegister;
    private int userDefaultQuota;

    private String databasePrefix;
    private String dataDir;
    private String tempJobsDir;
    private String toolsDir;

    private Project.Organism organism;

    private Admin admin;
    private List<AuthenticationOrigin> authenticationOrigins;
    private Monitor monitor;
    private Execution execution;
    private Audit audit;

    private List<StudyAclEntry> acl;

    private EmailServer emailServer;
    private DatabaseCredentials database;

    private RestServerConfiguration rest;
    private GrpcServerConfiguration grpc;

    protected static Logger logger = LoggerFactory.getLogger(Configuration.class);

    public Configuration() {
    }

//    public CatalogConfiguration(String defaultStorageEngineId, List<StorageEngineConfiguration> storageEngines) {
//        this.defaultStorageEngineId = defaultStorageEngineId;
//        this.storageEngines = storageEngines;
//
//        this.cellbase = new CellBaseConfiguration();
//        this.server = new QueryServerConfiguration();
//    }


    public void serialize(OutputStream configurationOututStream) throws IOException {
        ObjectMapper yamlMapper = new ObjectMapper(new YAMLFactory());
        yamlMapper.writerWithDefaultPrettyPrinter().writeValue(configurationOututStream, this);
    }

    public static Configuration load(InputStream configurationInputStream) throws IOException {
        return load(configurationInputStream, "yaml");
    }

    public static Configuration load(InputStream configurationInputStream, String format) throws IOException {
        Configuration configuration;
        ObjectMapper objectMapper;
        switch (format) {
            case "json":
                objectMapper = new ObjectMapper();
                configuration = objectMapper.readValue(configurationInputStream, Configuration.class);
                break;
            case "yml":
            case "yaml":
            default:
                objectMapper = new ObjectMapper(new YAMLFactory());
                configuration = objectMapper.readValue(configurationInputStream, Configuration.class);
                break;
        }
        return configuration;
    }

    @Override
    public String toString() {
        final StringBuilder sb = new StringBuilder("CatalogConfiguration{");
        sb.append("logLevel='").append(logLevel).append('\'');
        sb.append(", logFile='").append(logFile).append('\'');
        sb.append(", openRegister=").append(openRegister);
        sb.append(", userDefaultQuota=").append(userDefaultQuota);
        sb.append(", databasePrefix='").append(databasePrefix).append('\'');
        sb.append(", dataDir='").append(dataDir).append('\'');
        sb.append(", tempJobsDir='").append(tempJobsDir).append('\'');
        sb.append(", toolsDir='").append(toolsDir).append('\'');
        sb.append(", organism=").append(organism);
        sb.append(", admin=").append(admin);
        sb.append(", authenticationOrigins=").append(authenticationOrigins);
        sb.append(", monitor=").append(monitor);
        sb.append(", execution=").append(execution);
        sb.append(", audit=").append(audit);
        sb.append(", acl=").append(acl);
        sb.append(", emailServer=").append(emailServer);
        sb.append(", database=").append(database);
        sb.append('}');
        return sb.toString();
    }

    public String getLogLevel() {
        return logLevel;
    }

    public Configuration setLogLevel(String logLevel) {
        this.logLevel = logLevel;
        return this;
    }

    public String getLogFile() {
        return logFile;
    }

    public Configuration setLogFile(String logFile) {
        this.logFile = logFile;
        return this;
    }

    public boolean isOpenRegister() {
        return openRegister;
    }

    public Configuration setOpenRegister(boolean openRegister) {
        this.openRegister = openRegister;
        return this;
    }

    public int getUserDefaultQuota() {
        return userDefaultQuota;
    }

    public Configuration setUserDefaultQuota(int userDefaultQuota) {
        this.userDefaultQuota = userDefaultQuota;
        return this;
    }

    public String getDatabasePrefix() {
        return databasePrefix;
    }

    public Configuration setDatabasePrefix(String databasePrefix) {
        this.databasePrefix = databasePrefix;
        return this;
    }

    public String getDataDir() {
        return dataDir;
    }

    public Configuration setDataDir(String dataDir) {
        this.dataDir = dataDir;
        return this;
    }

    public String getTempJobsDir() {
        return tempJobsDir;
    }

    public Configuration setTempJobsDir(String tempJobsDir) {
        this.tempJobsDir = tempJobsDir;
        return this;
    }

    public String getToolsDir() {
        return toolsDir;
    }

    public Configuration setToolsDir(String toolsDir) {
        this.toolsDir = toolsDir;
        return this;
    }

    public Admin getAdmin() {
        return admin;
    }

    public Configuration setAdmin(Admin admin) {
        this.admin = admin;
        return this;
    }

    public List<AuthenticationOrigin> getAuthenticationOrigins() {
        return authenticationOrigins;
    }

    public Configuration setAuthenticationOrigins(List<AuthenticationOrigin> authenticationOrigins) {
        this.authenticationOrigins = authenticationOrigins;
        return this;
    }

    public Monitor getMonitor() {
        return monitor;
    }

    public Configuration setMonitor(Monitor monitor) {
        this.monitor = monitor;
        return this;
    }

    public Execution getExecution() {
        return execution;
    }

    public Configuration setExecution(Execution execution) {
        this.execution = execution;
        return this;
    }

    public EmailServer getEmailServer() {
        return emailServer;
    }

    public Configuration setEmailServer(EmailServer emailServer) {
        this.emailServer = emailServer;
        return this;
    }

    public DatabaseCredentials getDatabase() {
        return database;
    }

    public Configuration setDatabase(DatabaseCredentials database) {
        this.database = database;
        return this;
    }

    public Audit getAudit() {
        return audit;
    }

    public Configuration setAudit(Audit audit) {
        this.audit = audit;
        return this;
    }

    public List<StudyAclEntry> getAcl() {
        return acl;
    }

    public Configuration setAcl(List<StudyAclEntry> acl) {
        this.acl = acl;
        return this;
    }

    public Project.Organism getOrganism() {
        return organism;
    }

    public Configuration setOrganism(Project.Organism organism) {
        this.organism = organism;
        return this;
    }

    public RestServerConfiguration getRest() {
        return rest;
    }

    public Configuration setRest(RestServerConfiguration rest) {
        this.rest = rest;
        return this;
    }

    public GrpcServerConfiguration getGrpc() {
        return grpc;
    }

    public Configuration setGrpc(GrpcServerConfiguration grpc) {
        this.grpc = grpc;
        return this;
    }
}