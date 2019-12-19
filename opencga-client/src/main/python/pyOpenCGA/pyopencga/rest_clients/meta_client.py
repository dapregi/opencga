from pyopencga.rest_clients._parent_rest_clients import _ParentRestClient


class Meta(_ParentRestClient):
    """
    This class contains methods for the Meta webservice
    """

    def __init__(self, configuration, token=None, login_handler=None, *args, **kwargs):
        _category = "meta"
        super(Meta, self).__init__(configuration, _category, token, login_handler, *args, **kwargs)

    def about(self, **options):
        """
        Returns info about current OpenCGA code.
        URL: /{apiVersion}/meta/about
        """

        return self._get('about', **options)

    def ping(self, **options):
        """
        Ping Opencga webservices.
        URL: /{apiVersion}/meta/ping
        """

        return self._get('ping', **options)

    def status(self, **options):
        """
        Database status.
        URL: /{apiVersion}/meta/status
        """

        return self._get('status', **options)

    def api(self, **options):
        """
        API endpoints.
        URL: /{apiVersion}/meta/api
        """

        return self._get('api', **options)

    def fail(self, **options):
        """
        Ping Opencga webservices.
        URL: /{apiVersion}/meta/fail
        """

        return self._get('fail', **options)
