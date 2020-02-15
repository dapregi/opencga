import getpass
import time
import sys
import re

from pyopencga.opencga_config import ClientConfiguration
from pyopencga.rest_clients.admin_client import Admin
from pyopencga.rest_clients.alignment_client import Alignment
from pyopencga.rest_clients.clinical_client import Clinical
from pyopencga.rest_clients.cohort_client import Cohort
from pyopencga.rest_clients.family_client import Family
from pyopencga.rest_clients.file_client import File
from pyopencga.rest_clients.ga4gh_client import GA4GH
from pyopencga.rest_clients.individual_client import Individual
from pyopencga.rest_clients.job_client import Job
from pyopencga.rest_clients.meta_client import Meta
from pyopencga.rest_clients.disease_panel_client import DiseasePanel
from pyopencga.rest_clients.project_client import Project
from pyopencga.rest_clients.sample_client import Sample
from pyopencga.rest_clients.study_client import Study
from pyopencga.rest_clients.variant_operation_client import VariantOperation
from pyopencga.rest_clients.user_client import User
from pyopencga.rest_clients.variant_client import Variant


class OpencgaClient(object):
    def __init__(self, configuration, token=None, on_retry=None, auto_refresh=True):
        """
        :param on_retry: callback to be called with client retries an operation.
            It must accept parameters: client, exc_type, exc_val, exc_tb, call
        """

        if not isinstance(configuration, ClientConfiguration):
            raise ValueError('Expected ClientConfiguration instance')

        self.configuration = configuration
        self.auto_refresh = auto_refresh
        self.on_retry = on_retry
        self.clients = []
        self.user_id = None  # if user and session_id are supplied, we can log out
        self._login_handler = None
        self.token = token
        self._create_clients()
        self._check_versions()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.logout()

    def _check_versions(self):
        server_version = self.meta.about().get_result(0)['Version'].split('-')[0]
        client_version = re.findall(r'Client version: (.+)\n', str(self.meta.__doc__))[0]
        ansi_reset = "\033[0m"
        ansi_red = "\033[31m"
        ansi_yellow = "\033[33m"
        if tuple(server_version.split('.')[:2]) < tuple(client_version.split('.')[:2]):
            msg = '[WARNING]: Client version ({}) is higher than server version ({}).' \
                  ' Some client features may not be implemented in the server.\n'.format(client_version, server_version)
            sys.stdout.write('{}{}{}'.format(ansi_red, msg, ansi_reset))
        elif tuple(server_version.split('.')[:2]) > tuple(client_version.split('.')[:2]):
            msg = '[INFO]: Client version ({}) is lower than server version ({}).\n'.format(client_version, server_version)
            sys.stdout.write('{}{}{}'.format(ansi_yellow, msg, ansi_reset))

    def _create_clients(self):

        self.users = User(self.configuration, self.token, self._login_handler, auto_refresh=self.auto_refresh)
        self.projects = Project(self.configuration, self.token, self._login_handler, auto_refresh=self.auto_refresh)
        self.studies = Study(self.configuration, self.token, self._login_handler, auto_refresh=self.auto_refresh)
        self.files = File(self.configuration, self.token, self._login_handler, auto_refresh=self.auto_refresh)
        self.jobs = Job(self.configuration, self.token, self._login_handler, auto_refresh=self.auto_refresh)
        self.samples = Sample(self.configuration, self.token, self._login_handler, auto_refresh=self.auto_refresh)
        self.individuals = Individual(self.configuration, self.token, self._login_handler, auto_refresh=self.auto_refresh)
        self.families = Family(self.configuration, self.token, self._login_handler, auto_refresh=self.auto_refresh)
        self.cohorts = Cohort(self.configuration, self.token, self._login_handler, auto_refresh=self.auto_refresh)
        self.disease_panels = DiseasePanel(self.configuration, self.token, self._login_handler, auto_refresh=self.auto_refresh)
        self.alignments = Alignment(self.configuration, self.token, self._login_handler, auto_refresh=self.auto_refresh)
        self.variants = Variant(self.configuration, self.token, self._login_handler, auto_refresh=self.auto_refresh)
        self.clinical = Clinical(self.configuration, self.token, self._login_handler, auto_refresh=self.auto_refresh)
        self.variant_operations = VariantOperation(self.configuration, self.token, self._login_handler, auto_refresh=self.auto_refresh)
        self.meta = Meta(self.configuration, self.token, self._login_handler, auto_refresh=self.auto_refresh)
        self.ga4gh = GA4GH(self.configuration, self.token, self._login_handler, auto_refresh=self.auto_refresh)
        self.admin = Admin(self.configuration, self.token, self._login_handler, auto_refresh=self.auto_refresh)

        self.clients = [
            self.users, self.projects, self.studies, self.files, self.jobs,
            self.samples, self.individuals, self.families, self.cohorts,
            self.disease_panels, self.alignments, self.variants, self.clinical,
            self.variant_operations, self.meta, self.ga4gh, self.admin
        ]

        for client in self.clients:
            client.on_retry = self.on_retry

    def _make_login_handler(self, user, pwd):
        """
        Returns a closure that performs the log-in. This will be called on retries
        if the current session ever expires.
        The reason for using a closure and not a normal function is that a normal
        function would require storing the password in a field. It is more secure
        not to do so. This way, the password stored in the closure is inaccessible
        to other code
        """

        def login_handler(refresh=False):
            self.user_id = user
            if refresh:
                self.token = User(self.configuration).login(user=user, data={}).get_result(0)['token']
            else:
                self.token = User(self.configuration).login(user=user, data={'password': pwd}).get_result(0)['token']

            for client in self.clients:
                client.token = self.token  # renew the client's token
            return self.token

        return login_handler

    def login(self, user=None, password=None):
        if user is not None:
            if password is None:
                password = getpass.getpass()
            self._login_handler = self._make_login_handler(user, password)

        assert self._login_handler, "Can't login without username and password provided"
        self._login_handler()

    def logout(self):
        self.token = None
        for client in self.clients:
            client.token = self.token

    def wait_for_job(self, response=None, study=None, job_id=None, retry_seconds=10):
        if response is not None:
            study = response['studyUuid']
            job_id = response['uuid']

        if study is None or job_id is None:
            raise ValueError('Argument "response" or arguments "study" and "job_id" must be provided')

        if len(job_id.split(',')) > 1:
            raise ValueError('Only one job ID is allowed')

        retry_seconds = retry_seconds if retry_seconds >= 10 else 10
        while True:
            job_info = self.jobs.info(study=study, jobs=job_id, include='status').get_result(0)
            if job_info['status']['name'] in ['ERROR', 'ABORTED']:
                raise ValueError('{} ({}): {}'.format(
                    job_info['status']['name'], job_info['status']['date'], job_info['status']['message']
                ))
            elif job_info['status']['name'] in ['DONE']:
                break
            time.sleep(retry_seconds)

    def _get_help_info(self, method=None, parameters=False):
        info = []
        for client in self.clients:
            # Name
            class_name = type(client).__name__
            if class_name == 'GA4GH':
                class_name = class_name.lower()
            name = re.sub(r'(?<!^)(?=[A-Z])', '_', class_name).lower()
            name = 'get_' + name + '_client'

            if method is not None and method != name:
                continue

            # Description and path
            class_docstring = client.__doc__
            cls_desc = re.findall('(.+)\n +Client version', class_docstring)[0]
            cls_desc = cls_desc.strip().replace('This class contains methods',
                                                'Get client')
            cls_path = re.findall('PATH: (.+)\n', class_docstring)[0]

            # Methods
            methods = []
            method_names = [method_name for method_name in dir(client)
                            if callable(getattr(client, method_name))
                            and not method_name.startswith('_')]
            for method_name in method_names:
                method_docstring = getattr(client, method_name).__doc__
                desc = re.findall('(.+)\n +PATH', method_docstring, re.DOTALL)
                desc = re.sub(' +', ' ', desc[0].replace('\n', ' ').strip())
                path = re.findall('PATH: (.+)\n', method_docstring)[0]

                args = []
                arguments = re.findall(
                    ' +:param (.+)', method_docstring, re.DOTALL
                )
                if arguments and parameters:
                    arguments = arguments[0].replace('\n', ' ').strip()
                    arguments = re.sub(' +', ' ', arguments)
                    arguments = arguments.split(' :param ')
                    for parameter in arguments:
                        param_info = parameter.split(' ', 2)
                        args.append({
                            'name': param_info[1].rstrip(':'),
                            'type': param_info[0],
                            'desc': param_info[2]
                        })
                methods.append({
                    'name': method_name,
                    'desc': desc,
                    'path': path,
                    'params': args
                })

            info.append(
                {'name': name, 'desc': cls_desc, 'path': cls_path,
                 'methods': methods}
            )
        return info

    def help(self, method_name=None, show_parameters=False):
        help_text = []

        info = self._get_help_info(method_name, show_parameters)
        if method_name is None:
            for client_method in info:
                help_text += ['{}:'.format(client_method['name'])]
                help_text += ['{}- Description:'.format(' '*4)]
                help_text += ['{}{}'.format(' '*8, client_method['desc'])]
                help_text += ['{}- REST path:'.format(' '*4)]
                help_text += ['{}{}'.format(' '*8, client_method['path'])]
                help_text += ['{}- Usage:'.format(' '*4)]
                help_text += ['{}opencga_client.{}()\n'.format(' '*8, client_method['name'])]
        else:
            for client_method in info:
                if method_name == client_method['name']:
                    help_text.append(client_method['name'] + '\n')
                    help_text += [m['name'] for method in info for m in method['methods']]

        sys.stdout.write('\n'.join(help_text) + '\n')


    def get_user_client(self):
        return self.users

    def get_project_client(self):
        return self.projects

    def get_study_client(self):
        return self.studies

    def get_file_client(self):
        return self.files

    def get_job_client(self):
        return self.jobs

    def get_sample_client(self):
        return self.samples

    def get_individual_client(self):
        return self.individuals

    def get_family_client(self):
        return self.families

    def get_cohort_client(self):
        return self.cohorts

    def get_disease_panel_client(self):
        return self.disease_panels

    def get_alignment_client(self):
        return self.alignments

    def get_variant_client(self):
        return self.variants

    def get_clinical_client(self):
        return self.clinical

    def get_variant_operation_client(self):
        return self.variant_operations

    def get_meta_client(self):
        return self.meta

    def get_ga4gh_client(self):
        return self.ga4gh

    def get_admin_client(self):
        return self.admin
