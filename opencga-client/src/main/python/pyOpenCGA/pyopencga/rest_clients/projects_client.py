from pyopencga.rest_clients._parent_rest_clients import _ParentRestClient


class Projects(_ParentRestClient):
    """
    This class contains methods for the 'Projects' webservices
    Client version: 2.0.0
    PATH: /{apiVersion}/projects
    """

    def __init__(self, configuration, token=None, login_handler=None, *args, **kwargs):
        _category = 'projects'
        super(Projects, self).__init__(configuration, _category, token, login_handler, *args, **kwargs)

    def update(self, project, data, **options):
        """
        Update some project attributes
        PATH: /{apiVersion}/projects/{project}/update

        :param str project: Project [user@]project where project can be either the ID or the alias
        :param dict data: JSON containing the params to be updated. It will be only possible to update organism fields not previously defined.
        """

        return self._post('update', query_id=project, data=data, **options)

    def aggregation_stats(self, projects, **options):
        """
        Fetch catalog project stats
        PATH: /{apiVersion}/projects/{projects}/aggregationStats

        :param str projects: Comma separated list of projects [user@]project up to a maximum of 100
        :param bool default: Calculate default stats
        :param str file_fields: List of file fields separated by semicolons, e.g.: studies;type. For nested fields use >>, e.g.: studies>>biotype;type
        :param str individual_fields: List of individual fields separated by semicolons, e.g.: studies;type. For nested fields use >>, e.g.: studies>>biotype;type
        :param str family_fields: List of family fields separated by semicolons, e.g.: studies;type. For nested fields use >>, e.g.: studies>>biotype;type
        :param str sample_fields: List of sample fields separated by semicolons, e.g.: studies;type. For nested fields use >>, e.g.: studies>>biotype;type
        :param str cohort_fields: List of cohort fields separated by semicolons, e.g.: studies;type. For nested fields use >>, e.g.: studies>>biotype;type
        """

        return self._get('aggregationStats', query_id=projects, **options)

    def create(self, data, **options):
        """
        Create a new project
        PATH: /{apiVersion}/projects/create

        :param dict data: JSON containing the mandatory parameters
        """

        return self._post('create', data=data, **options)

    def search(self, **options):
        """
        Search projects
        PATH: /{apiVersion}/projects/search

        :param str include: Fields included in the response, whole JSON path must be provided
        :param str exclude: Fields excluded in the response, whole JSON path must be provided
        :param int limit: Number of results to be returned
        :param int skip: Number of results to skip
        :param str owner: Owner of the project
        :param str id: Project [user@]project where project can be either the ID or the alias
        :param str name: Project name
        :param str fqn: Project fqn
        :param str alias: DEPRECATED: Project alias
        :param str organization: Project organization
        :param str description: Project description
        :param str study: Study id or alias
        :param str creation_date: Creation date. Format: yyyyMMddHHmmss. Examples: >2018, 2017-2018, <201805
        :param str modification_date: Modification date. Format: yyyyMMddHHmmss. Examples: >2018, 2017-2018, <201805
        :param str status: Status
        :param str attributes: Attributes
        """

        return self._get('search', **options)

    def inc_release(self, project, **options):
        """
        Increment current release number in the project
        PATH: /{apiVersion}/projects/{project}/incRelease

        :param str project: Project [user@]project where project can be either the ID or the alias
        """

        return self._post('incRelease', query_id=project, **options)

    def studies(self, project, **options):
        """
        Fetch all the studies contained in the project
        PATH: /{apiVersion}/projects/{project}/studies

        :param str include: Fields included in the response, whole JSON path must be provided
        :param str exclude: Fields excluded in the response, whole JSON path must be provided
        :param int limit: Number of results to be returned
        :param int skip: Number of results to skip
        :param str project: Project [user@]project where project can be either the ID or the alias
        """

        return self._get('studies', query_id=project, **options)

    def info(self, projects, **options):
        """
        Fetch project information
        PATH: /{apiVersion}/projects/{projects}/info

        :param str include: Fields included in the response, whole JSON path must be provided
        :param str exclude: Fields excluded in the response, whole JSON path must be provided
        :param str projects: Comma separated list of projects [user@]project up to a maximum of 100
        """

        return self._get('info', query_id=projects, **options)
