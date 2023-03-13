# =========================================================
# ImExULTIMO.py
# @author Nina Thomsen
# @date 28.09.2022
# @copyright Institut fuer Verkehrsforschung,
#            Deutsches Zentrum fuer Luft- und Raumfahrt
# @brief Sets type of and executes import, export of files
# =========================================================

import psycopg2
import geopandas as gpd
import pandas as pd
from sqlalchemy import create_engine


class ImportU:

    """This class is used for importing data into ULTIMO"""

    def __init__(self, _type='local'):
        self.type = _type

    def read_geo_input(self, user="", password="", host="", port=0, dbname="", query="", filename="", geom="geom"):
        """
        :param port: Type postgres - port
        :param dbname: Type postgres - database name
        :param host: Type postgres - database host
        :param password: Type postgres - password
        :param user: Type postgres - username
        :param filename: Type local - location of file in directory
        :param query: Type postgres - query for data from postgres
        :param geom: Type postgres - name of geometry column
        """
        input_ = False
        if self.type == "postgres":
            engine = create_engine('postgresql+psycopg2://{}:{}@{}:{}/{}'.format(user, password, host, port, dbname))
            input_ = gpd.GeoDataFrame.from_postgis(query, engine, geom_col=geom)

        if self.type == "local":
            input_ = gpd.GeoDataFrame.from_file(filename)

        return input_

    def read_table_input(self, postgres="", query="", filename="", sep=","):
        """
        :param sep: Type local - separator of CSV file
        :param postgres: Type postgres - connection data for database as str
        :param filename: Type local - location of file in directory
        :param query: Type postgres - query for data from postgres
        """
        input_ = False
        if self.type == "postgres":
            conn = psycopg2.connect(postgres)
            input_ = pd.read_sql(query, conn)

        if self.type == "local":
            input_ = pd.read_csv(filename, sep=sep)

        return input_


class ExportU:

    """
    This class is used for exporting data from any ULTIMO module.
    Use with:
     type 'postgres' for database
     type 'local' for local file
     default 'local'
    """

    def __init__(self, _file, _type='local'):
        self.type = _type
        self.file = _file

    def geo_output(self, user="", password="", host="", port=0, dbname="", filename="", schema="", path_=""):
        """
        :param schema: Type postgres - schema name
        :param path_: Type local - directory for output
        :param port: Type postgres - port
        :param dbname: Type postgres - database name
        :param host: Type postgres - database host
        :param password: Type postgres - password
        :param user: Type postgres - username
        :param filename: Name of table / file
        """

        if self.type == "postgres":
            engine = create_engine('postgresql+psycopg2://{}:{}@{}:{}/{}'.format(user, password, host, port, dbname))
            self.file.to_postgis(filename, engine, schema=schema, if_exists='replace')

        if self.type == "local":
            self.file.to_file(path_+filename+".gpkg", driver="GPKG")

    def table_output(self, user="", password="", host="", port=0, dbname="", filename="", schema="", path_=""):
        """
        :param schema: Type postgres - schema name
        :param path_: Type local - directory for output
        :param port: Type postgres - port
        :param dbname: Type postgres - database name
        :param host: Type postgres - database host
        :param password: Type postgres - password
        :param user: Type postgres - username
        :param filename: Name of table / file
        """

        if self.type == "postgres":
            engine = create_engine('postgresql+psycopg2://{}:{}@{}:{}/{}'.format(user, password, host, port, dbname))
            self.file.to_sql(filename, engine, schema=schema, if_exists='replace')

        if self.type == "local":
            self.file.to_csv(path_+filename+".csv")
