from pyopencga.rest_clients._parent_rest_clients import _ParentRestClient


class Panels(_ParentRestClient):
    """
    This class contains methods for the 'Disease Panels' webservices
    Client version: 2.0.0
    PATH: /{apiVersion}/panels
    """

    def __init__(self, configuration, token=None, login_handler=None, *args, **kwargs):
        _category = 'panels'
        super(Panels, self).__init__(configuration, _category, token, login_handler, *args, **kwargs)

    def update(self, panels, data=None, **options):
        """
        Update panel attributes
        PATH: /{apiVersion}/panels/{panels}/update

        :param str study: Study [[user@]project:]study where study and project can be either the ID or UUID
        :param str panels: Comma separated list of panel ids
        :param bool inc_version: Create a new version of panel
        :param dict data: Panel parameters
        """

        return self._post('update', query_id=panels, data=data, **options)

    def delete(self, panels, **options):
        """
        Delete existing panels
        PATH: /{apiVersion}/panels/{panels}/delete

        :param str study: Study [[user@]project:]study where study and project can be either the ID or UUID
        :param str panels: Comma separated list of panel ids
        """

        return self._delete('delete', query_id=panels, **options)

    def create(self, data=None, **options):
        """
        Create a panel
        PATH: /{apiVersion}/panels/create

        :param str study: Study [[user@]project:]study where study and project can be either the ID or UUID
        :param str import: Comma separated list of installation panel ids to be imported. To import them all at once, write the special word 'ALL_GLOBAL_PANELS'
        :param dict data: Panel parameters
        """

        return self._post('create', data=data, **options)

    def acl(self, panels, **options):
        """
        Returns the acl of the panels. If member is provided, it will only return the acl for the member.
        PATH: /{apiVersion}/panels/{panels}/acl

        :param str panels: Comma separated list of panel ids up to a maximum of 100
        :param str study: Study [[user@]project:]study where study and project can be either the ID or UUID
        :param str member: User or group id
        :param bool silent: Boolean to retrieve all possible entries that are queried for, false to raise an exception whenever one of the entries looked for cannot be shown for whichever reason
        """

        return self._get('acl', query_id=panels, **options)

    def update_acl(self, members, data, **options):
        """
        Update the set of permissions granted for the member
        PATH: /{apiVersion}/panels/acl/{members}/update

        :param str study: Study [[user@]project:]study where study and project can be either the ID or UUID
        :param str members: Comma separated list of user or group ids
        :param dict data: JSON containing the parameters to update the permissions.
        """

        return self._post('update', query_id=members, data=data, **options)

    def info(self, panels, **options):
        """
        Panel info
        PATH: /{apiVersion}/panels/{panels}/info

        :param str include: Fields included in the response, whole JSON path must be provided
        :param str exclude: Fields excluded in the response, whole JSON path must be provided
        :param str panels: Comma separated list of panel ids up to a maximum of 100
        :param str study: Study [[user@]project:]study where study and project can be either the ID or UUID
        :param int version: Panel  version
        :param bool deleted: Boolean to retrieve deleted panels
        :param bool global: Boolean indicating which panels are queried (installation or study specific panels)
        """

        return self._get('info', query_id=panels, **options)

    def search(self, **options):
        """
        Panel search
        PATH: /{apiVersion}/panels/search

        :param str include: Fields included in the response, whole JSON path must be provided
        :param str exclude: Fields excluded in the response, whole JSON path must be provided
        :param int limit: Number of results to be returned
        :param int skip: Number of results to skip
        :param bool count: Get the total number of results matching the query. Deactivated by default.
        :param str study: Study [[user@]project:]study where study and project can be either the ID or UUID
        :param str name: Panel name
        :param str phenotypes: Panel phenotypes
        :param str variants: Panel variants
        :param str genes: Panel genes
        :param str regions: Panel regions
        :param str categories: Panel categories
        :param str tags: Panel tags
        :param str description: Panel description
        :param str author: Panel author
        :param bool deleted: Boolean to retrieve deleted panels
        :param str creation_date: Creation date. Format: yyyyMMddHHmmss. Examples: >2018, 2017-2018, <201805
        :param str modification_date: Modification date. Format: yyyyMMddHHmmss. Examples: >2018, 2017-2018, <201805
        :param bool global: Boolean indicating which panels are queried (installation or study specific panels)
        :param str release: Release value (Current release from the moment the samples were first created)
        :param int snapshot: Snapshot value (Latest version of samples in the specified release)
        """

        return self._get('search', **options)
