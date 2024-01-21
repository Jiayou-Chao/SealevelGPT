import openai
from openai import OpenAI
from loguru import logger
import os
import sys
import utils
import subprocess


def set_log_level(level=None):
    """Set the log level for the logger."""
    if level is None:
        level = os.getenv('LOG_LEVEL', 'DEBUG')
    logger.remove()
    logger.add(
        sys.stderr,
        level=level,
    )
    

def find_location_and_year_with_gpt(question: list | str = None, 
                                    model="gpt-3.5-turbo", 
                                    # chat_id=None, 
                                    conversation_messages=None,
                                    max_trials=1,
                                    nth_trial=0,
                                    **kwargs):

    logger.debug(f"Model: {model}")
    logger.info(f"User messages: {question}")

    if conversation_messages is None:
        system_message = (
            "The assistant is a highly capable entity detector. "
            "It identifies and extracts locations and years from user inputs. "
            "The user inputs should be in the form of a question that asks for the sea level prediction for a location and year."
            "For instance, if the user asks 'what is the population of Tokyo in 2021?', "
            "the assistant should respond 'Location: Tokyo, Year: 2021'."
            "If you think the question is not clear enough or is not related to sea level, you should express your concern and ask the user to rephrase the question."
            "Your response should be as succinct as possible."
        )
        conversation_messages = [
            {"role": "system", "content": system_message}
        ]

    if question is None:
        question = input("Please enter a question that asks for the sea level prediction for a location and year: ")
    if not isinstance(question, list):
        question = [question]
    
    # Add user messages to the conversation
    conversation_messages.extend(
        {"role": "user", "content": msg} for msg in question
    )

    # # Append an instruction to analyze the messages and summarize
    # conversation_messages.append({
    #     "role": "system",
    #     "content": "Please read the above conversation and identify any mentioned locations and years."
    # })

    # Create the Chat object and send the messages
    client = OpenAI(
        # This is the default and can be omitted
        api_key=utils.get_api_key(),
    )
    response = client.chat.completions.create(
        messages=conversation_messages,
        model=model,
        **kwargs
    )

    # Extract the assistant's message
    assistant_message = response.choices[0].message.content
    # new_chat_id = response['data']['chat']
    # logger.debug(f"New Chat ID: {new_chat_id}")
    logger.debug(f"conversation_messages: {conversation_messages}")
    logger.debug(f"Assistant message: {assistant_message}")

    # Check if the response contains information about location and year
    if "Location:" not in assistant_message and "Year:" not in assistant_message:
        nth_trial += 1
        if nth_trial < max_trials:
            logger.info("Re-trying to find the location and year.")
            return find_location_and_year_with_gpt(question=None, model=model, conversation_messages=conversation_messages, max_trials=max_trials, nth_trial=nth_trial, **kwargs)
        else:
            raise ValueError(f"Cannot find location and year after {max_trials} trials.")

    location = assistant_message.split("Location:")[1].split(",")[0].strip()
    year = assistant_message.split("Year:")[1].strip()
    logger.info(f"Location: {location}, Year: {year}")
    return location, year


def run_r_script(location, year):
    R_script = f"""
    source("chatgpt-sealevel.R")
    GPT_predict("{location}", "{year}")
    """
    with open("temporary_chatgpt-sealevel.R", "w") as f:
        f.write(R_script)
    # os.system("Rscript temporary_chatgpt-sealevel.R > /dev/null 2>&1")
    with open(os.devnull, 'w') as devnull:
        subprocess.run(["Rscript", "temporary_chatgpt-sealevel.R"], stdout=devnull, stderr=subprocess.STDOUT)
    with open("temporary_prediction.txt", "r") as f:
        output = f.read()
    os.remove("temporary_prediction.txt")
    os.remove("temporary_chatgpt-sealevel.R")
    return output
    

def predict_local_sea_level_with_gpt(question = None, model="gpt-3.5-turbo", **kwargs):
    if question is None:
        question = input("Please enter a question that asks for the sea level prediction for a location and year: ")
    location, year = find_location_and_year_with_gpt(question, model=model)
    location = utils.format_location(location)
    
    output = run_r_script(location, year)
    logger.debug(f"Output: {output}")

    system_message = (
        "Answer the user's question based on the following information: " +
        output +
        " You should also give analysis related to global warming and sea level rise based on the prediction."
    )

    conversation_messages = [
        {"role": "system", "content": system_message},
        {"role": "user", "content": question}
    ]

    client = OpenAI(
        # This is the default and can be omitted
        api_key=utils.get_api_key(),
    )
    response = client.chat.completions.create(
        messages=conversation_messages,
        model="gpt-4-1106-preview",
        # model="gpt-3.5-turbo",
        **kwargs
    )

    assistant_message = response.choices[0].message.content
    print(assistant_message)
    conversation_messages.append(
        {"role": "assistant", "content": assistant_message}
    )
    i = 0
    while question != 'stop' and question != 'quit' and i < 5:
        i += 1
        question = input("Do you have any other questions? Type 'stop' or 'quit' or press Ctrl+C to end the conversation: ")
        conversation_messages.append(
            {"role": "user", "content": question}
        )
        response = client.chat.completions.create(
            messages=conversation_messages,
            model="gpt-4-1106-preview",
            **kwargs
        )
        conversation_messages.append(
            {"role": "assistant", "content": response.choices[0].message.content}
        )
    logger.debug(f"conversation_messages: {conversation_messages}")
    
    with open("chatgpt-response.txt", "w") as f:
        f.write(assistant_message)
    return assistant_message
        
    


def main():
    # Set the log level
    set_log_level('DEBUG')
    # Test the function
    # question = "What is the sea level in tofino in 2050?"
    # What is the sea level in tofino in 2050?
    # predict_local_sea_level_with_gpt(question, model="gpt-3.5-turbo")
    predict_local_sea_level_with_gpt(model="gpt-3.5-turbo")
    
if __name__ == "__main__":
    main()
        
        
    
    
    
