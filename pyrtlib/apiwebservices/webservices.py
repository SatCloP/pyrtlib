"""Utility code to support making requests using HTTP.
The orignal code of this class can be found in https://github.com/Unidata/siphon
Terms and conditions are as for siphon library license (https://github.com/Unidata/siphon/blob/master/LICENSE)
"""

import posixpath
from io import BytesIO
from typing import Optional

import requests
from requests.sessions import Session

from ..version import __version__

HTTPError = requests.HTTPError


class HTTPSessionManager(object):
    """Manage the creation of sessions for HTTP access."""

    def __init__(self):
        """Initialize ``HTTPSessionManager``."""
        self.user_agent = f'pyrtlib ({__version__})'
        self.options = {}

    def set_session_options(self, **kwargs) -> None:
        """Set options for created session instances.

        Takes keyword arguments and sets them as attributes on the returned
        :class:`requests.Session` instance.

        See Also
        --------
        create_session

        """
        self.options = kwargs

    def create_session(self) -> Session:
        """Create a new HTTP session with our user-agent set.

        Returns
        -------
        session : requests.Session
            The created session

        See Also
        --------
        urlopen, set_session_options

        """
        ret = requests.Session()
        ret.headers['User-Agent'] = self.user_agent
        for k, v in self.options.items():
            setattr(ret, k, v)
        return ret

    def urlopen(self, url: str, **kwargs) -> BytesIO:
        """GET a file-like object for a URL using HTTP.

        This is a thin wrapper around :meth:`requests.Session.get` that returns a file-like
        object wrapped around the resulting content.

        Parameters
        ----------
        url : str
            The URL to request

        kwargs : arbitrary keyword arguments
            Additional keyword arguments to pass to :meth:`requests.Session.get`.

        Returns
        -------
        fobj : file-like object
            A file-like interface to the content in the response

        See Also
        --------
        :meth:`requests.Session.get`

        """
        return BytesIO(self.create_session().get(url, **kwargs).content)


session_manager = HTTPSessionManager()


class HTTPEndPoint(object):
    """An object representing an endpoint on a server that is accessed using HTTP.

    This provides a simple way to point to a URL, formulate appropriate queries and
    validate them, parse metadata as appropriate, and parse returns from requests.
    """

    def __init__(self, url: str) -> None:
        """Create an HTTPEndPoint instance.

        Parameters
        ----------
        url : str
            The base URL for the endpoint

        """
        self._base = url
        self._session = session_manager.create_session()
        self._get_metadata()

    def _get_query(self, query: dict) -> str:
        """Make a GET request, including a query, to the endpoint.

        The path of the request is to the base URL assigned to the endpoint.

        Parameters
        ----------
        query : dict
            The query to pass when making the request

        Returns
        -------
        resp : requests.Response
            The server's response to the request

        See Also
        --------
        get_path, get

        """
        url = self._base[:-1] if self._base[-1] == '/' else self._base
        return self._get(url, query)

    def _url_path(self, path: str) -> str:
        """Assemble the full url to a path.

        Given a path relative to the base URL, assemble the full URL.

        Parameters
        ----------
        path : str
            The path, relative to the endpoint

        Returns
        -------
        url : str
            The full URL to `path`

        See Also
        --------
        get_path

        """
        return posixpath.join(self._base, path)

    def _get_path(self, path: str, querystrings: Optional[dict] = None) -> requests.Response:
        """Make a GET request, optionally including a query, to a relative path.

        The path of the request includes a path on top of the base URL
        assigned to the endpoint.

        Parameters
        ----------
        path : str
            The path to request, relative to the endpoint
        querystrings : dict, optional
            The query to pass when making the request

        Returns
        -------
        resp : requests.Response
            The server's response to the request

        See Also
        --------
        get_query, get, url_path

        """
        return self._get(self._url_path(path), querystrings)

    def _get(self, path: str, params: Optional[dict] = None) -> str:
        """Make a GET request, optionally including a parameters, to a path.

        The path of the request is the full URL.

        Parameters
        ----------
        path : str
            The URL to request
        params : dict, optional
            The query to pass when making the request

        Returns
        -------
        resp : requests.Response
            The server's response to the request

        Raises
        ------
        HTTPError
            If the server returns anything other than a 200 (OK) code

        See Also
        --------
        get_query, get

        """
        resp = self._session.get(path, params=params)
        if resp.status_code != 200:
            if resp.headers.get('Content-Type', '').startswith('text/html'):
                text = resp.reason
            else:
                text = resp.text
            raise requests.HTTPError('Error accessing {0}\n'
                                     'Server Error ({1:d}: {2})'.format(resp.request.url,
                                                                        resp.status_code,
                                                                        text))
        return resp

    def _get_metadata(self):
        """Get the metadata associated with the endpoint.

        It is intended that this be implemented by subclasses as necessary.
        """
