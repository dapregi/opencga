package org.opencb.opencga.catalog.managers.api;

import org.opencb.commons.datastore.core.ObjectMap;
import org.opencb.commons.datastore.core.QueryOptions;
import org.opencb.commons.datastore.core.QueryResult;
import org.opencb.opencga.catalog.exceptions.CatalogException;
import org.opencb.opencga.catalog.models.QueryFilter;
import org.opencb.opencga.catalog.models.User;

import java.io.IOException;

/**
 * @author Jacobo Coll &lt;jacobo167@gmail.com&gt;
 */
public interface IUserManager extends ResourceManager<String, User> {

    /**
     * Get the userId from the sessionId.
     *
     * @param sessionId SessionId
     * @return UserId owner of the sessionId. Empty string if SessionId does not match.
     * @throws CatalogException when the session id does not correspond to any user or the token has expired.
     */
    String getUserId(String sessionId) throws CatalogException;

    /**
     * Create a new user.
     *
     * @param userId       User id
     * @param name         Name
     * @param email        Email
     * @param password     Encrypted Password
     * @param organization Optional organization
     * @param diskQuota    Maximum user disk quota
     * @param options      Optional options
     * @param adminPassword    Optional admin password.
     * @return The created user
     * @throws CatalogException If user already exists, or unable to create a new user.
     */
    QueryResult<User> create(String userId, String name, String email, String password, String organization, Long diskQuota,
                             QueryOptions options, String adminPassword) throws CatalogException;

    /**
     * Gets the user information.
     *
     * @param userId       User id
     * @param lastModified If lastModified matches with the one in Catalog, return an empty QueryResult.
     * @param options      QueryOptions
     * @param sessionId    SessionId of the user performing this operation.
     * @return The requested user
     * @throws CatalogException CatalogException
     */
    QueryResult<User> read(String userId, String lastModified, QueryOptions options, String sessionId) throws CatalogException;

    void changePassword(String userId, String oldPassword, String newPassword) throws CatalogException;

    @Deprecated
    QueryResult<ObjectMap> login(String userId, String password, String sessionIp) throws CatalogException, IOException;

    /**
     * This method will be only callable by the admin. It generates a new session id for the user.
     *
     * @param sessionId Admin session id.
     * @param userId user id for which a session will be generated.
     * @return an objectMap containing the new sessionId
     * @throws CatalogException if the password is not correct or the userId does not exist.
     */
    QueryResult<ObjectMap> getNewUserSession(String sessionId, String userId) throws CatalogException;

    QueryResult resetPassword(String userId) throws CatalogException;

    @Deprecated
    QueryResult<ObjectMap> loginAsAnonymous(String sessionIp) throws CatalogException, IOException;

    QueryResult logout(String userId, String sessionId) throws CatalogException;

    @Deprecated
    QueryResult logoutAnonymous(String sessionId) throws CatalogException;

    void addQueryFilter(String sessionId, QueryFilter queryFilter) throws CatalogException;

    QueryResult<Long> deleteQueryFilter(String sessionId, String filterId) throws CatalogException;

    QueryResult<QueryFilter> getQueryFilter(String sessionId, String filterId) throws CatalogException;

}
