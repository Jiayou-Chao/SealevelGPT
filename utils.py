import openai
from openai import OpenAI
from dotenv import load_dotenv
import os
import re

def format_location(location):
    """
    remove leading and trailing spaces and convert to lowercase. Replace the spaces with underscores.
    For example, ' Astoria ' becomes 'astoria'.
    """
    location = location.strip().lower().replace(" ", "_")
    return location

def format_year(year: str):
    """
    remove any non-digit characters from the year and convert to an integer.
    For example, ' 2021.' becomes 2021.
    """
    year = re.sub("[^0-9]", "", year)
    year = int(year)
    return year

def get_api_key():
    """
    Get the API key from the .env file.
    """
    load_dotenv()
    api_key = os.getenv("OPENAI_API_KEY")
    return api_key

def list_available_models():
    """
    List all the models available on OpenAI.
    """
    client = OpenAI(
        api_key=get_api_key()
    )
    models = client.models.list()
    return models

if __name__ == "__main__":
    print(list_available_models())



    