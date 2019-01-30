#!/usr/bin/python

import re
import sys


def make_class_header(class_name, description, space=""):

    # remove // from description
    description = re.sub(r"//", "", description).strip()
    description = re.sub(r"\n[ ]*", "\n* ", description)
    class_name = re.sub(r"{", "", class_name).strip()
    class_name = class_name.split(':')[0].strip()

    boilerplate = f"""{space}
{space}/**
{space}* \\brief {description}
{space}*/"""

    return boilerplate


def make_method_header(description="", parameters=[], space=""):
    # remove // from description
    description = re.sub(r"//", "", description).strip()
    description = re.sub(r"\n[ ]*", "\n* ", description)

    boilerplate = ""

    if description != "":

        boilerplate += f"""{space}
{space}/**
{space}* \\brief {description}
{space}*
"""
    elif parameters != []:
        boilerplate += f"""{space}
{space}/**
"""

    if parameters != []:
        for param in parameters:
            boilerplate += f"""{space}* \param {(param.split('=')[0].strip()).split(' ')[-1]}
"""

        boilerplate += f"""{space}*/
"""

    return boilerplate


def make_method_doxycomment(description="", space=""):
    # remove // from description
    description = re.sub(r"//", "", description).strip()
    description = re.sub(r"\n[ ]*", "\n\* ", description)

    if description == "":
        return ""

    else:

        return f"""{space}
{space}/**
{space}* \\note
{space}* {description}
{space}*/
"""


def make_variable_docstring(description, space=""):
    description = re.sub(r"//", "", description).strip()
    description = re.sub(r"\n[ ]*", "\n* ", description)

    if description == "":
        return ""
    else:
        return f"""{space}
{space}/**
{space}* {description}
{space}*/
"""


def process_header_file(filename):

    output_data = ""

    # find comments in lines above
    re_comments = re.compile(
        r"([ \t]*)\/\/\s*\n((?:[ \t]*\/\/[\S ]*\n)+)[ \t]*\/\/\s*$", flags=re.MULTILINE)

    with open(filename) as input_file:
        data = input_file.read()

        # find classes
        re_class_name = re.compile(r"\n[^\n\s]*class ([^\n]+)")

        last_index = 0

        for m in re.finditer(re_class_name, data):

            comments = None
            for comments in re.finditer(re_comments,
                                        data[last_index:m.start()]):
                pass

            if comments and (m.start() - comments.end() - last_index) < 3:
                output_data += data[last_index:last_index + comments.start()]
                class_header = make_class_header(m.group(1), comments.group(2), space=comments.group(1))
            else:
                output_data += data[last_index:m.start()]
                class_header = make_class_header(m.group(1), "")

            output_data += class_header

            last_index = m.start()

        output_data += data[last_index:]

    data = output_data
    output_data = ""
    last_index = 0

    re_prototype = re.compile(
        r"(?:^[\w&:*\t ]+\n)*^([ \t]*)[~\w:*& <>]+\(([*\w\: \,&\n\t_=\<>\-.]*)\)", flags=re.MULTILINE)

    # markup methods
    for m in re.finditer(re_prototype, data):
        # print("match = ", m.group(1))

        parameters = m.group(2).split(",")
        parameters = [param.strip() for param in parameters]
        parameters = [param for param in parameters if param != ""]

        comments = None
        for comments in re.finditer(re_comments,
                                    data[last_index:m.start()]):
            pass

        if comments and (m.start() - comments.end() - last_index) < 2:
            # print(comments.span())
            output_data += data[last_index:last_index + comments.start()]
            method_header = make_method_header(comments.group(2), parameters, space=m.group(1))
            last_index = m.start()

        else:
            output_data += data[last_index:m.start()]
            method_header = make_method_header("", parameters, space=m.group(1))
            last_index = m.start()

        output_data += method_header

    output_data += data[last_index:]

    data = output_data
    output_data = ""
    last_index = 0

    re_comments = re.compile(
        r"^[ \t]*(\/\/[ \t]*[\S \n]*?)\n^(?![ \t]*\/\/)", flags=re.MULTILINE)

    re_variable = re.compile(
        r"^([ \t]*)[~\w:*& <>\[\]]+;", flags=re.MULTILINE)

    # markup variables
    for m in re.finditer(re_variable, data):

        if " return " in m.group(0):
            continue

        comments = None
        for comments in re.finditer(re_comments,
                                    data[last_index:m.start() + 1]):
            pass

        if comments and (m.start() - comments.end() - last_index) < 1:
            output_data += data[last_index:last_index + comments.start()]
            variable_header = make_variable_docstring(comments.group(1), space=m.group(1))
            output_data += variable_header
            last_index = m.start()

        else:
            output_data += data[last_index:m.start()]
            last_index = m.start()

    output_data += data[last_index:]

    output_filename = filename + ".doxygen"

    with open(output_filename, 'w+') as output_file:
        output_file.write(output_data)


def process_cpp_file(filename):

    output_data = ""

    # find comments in lines above
    re_comments = re.compile(
        r"[ \t]*(\/\/[ \t]*[\S ^\n]*?)\n^[ \t]*[^\/]", flags=re.MULTILINE)

    re_prototype = re.compile(
        r"^\w*\n^\w[~\w:*& ]+\([\w\: \,&\n\t_=<>.]*\)\n?[\s\S]*?\n?{", flags=re.MULTILINE)

    with open(filename) as input_file:
        data = input_file.read()

        last_index = 0

        for m in re.finditer(re_prototype, data):

            comments = None

            for comments in re.finditer(re_comments,
                                        data[last_index:m.start() + 1]):
                pass

            if comments and (m.start() - comments.end() - last_index) < 3:
                output_data += data[last_index:last_index + comments.start()]
                method_header = make_method_doxycomment(comments.group(1))
                last_index = m.start()

            else:
                output_data += data[last_index:m.start()]
                method_header = make_method_doxycomment("")
                last_index = m.start()

            output_data += method_header

        output_data += data[last_index:]

    output_filename = filename + ".doxygen"

    with open(output_filename, 'w+') as output_file:
        output_file.write(output_data)

def make_subroutine_header(description="", binds_to="", parameters=[]):
    description = re.sub(r"!", "", description).strip()
    description = re.sub(r"\n[ ]*", "\n!! ", description)

    if description == "" and binds_to == "" and parameters == []:
        return ""

    boilerplate = """
!> """
    if description != "":
        boilerplate += f"""\brief {description}
!!
"""
    else:
        boilerplate += """
"""

    if binds_to != "":
        boilerplate += f"""!! \note Binds to C function ``{binds_to.strip()}``
!!
"""

    if parameters != []:

        for (type, intent, vars) in parameters:
            vars = vars.split(",")
            for v in vars:
                boilerplate += f"""!! \param[{intent.strip()}] {v.strip()} {type.strip()}
"""

        boilerplate += """!!
"""

    return boilerplate

def make_function_header(description="", parameters=[]):
    description = re.sub(r"!", "", description).strip()
    description = re.sub(r"\n[ ]*", "\n!! ", description)

    if description == "" and parameters == []:
        return ""

    boilerplate = """
!> """
    if description != "":
        boilerplate += f"""\brief {description}
!!
"""
    else:
        boilerplate += """
"""
    if parameters != []:

        for (type, intent, vars) in parameters:
            vars = vars.split(",")
            for v in vars:
                boilerplate += f"""!! \param[{intent.strip()}] {v.strip()} {type.strip()}
"""

        boilerplate += """!!
"""

    return boilerplate


def process_fortran_file(filename):
    output_data = ""

    re_subroutine = re.compile(r"^[ \t]*(?:pure )?subroutine \S+\([\w ,&\n]*\)(?:[ &\n]*bind\(C, *name *= ?\"(\w+)\")?", flags=re.MULTILINE)

    re_end_subroutine = re.compile(r"^[ \t]*end subroutine", flags=re.MULTILINE)

    re_comments = re.compile(r"^[ \t]*!([\S \n]*?)(?=^[ \t]*[^!]*$)", flags=re.MULTILINE)

    re_parameters = re.compile(r"^[ \t]*([\w() ]+) *, *intent *\( *(in|out|inout) *\) *:: *([\w, ]+)", flags=re.MULTILINE)

    with open(filename) as input_file:
        data = input_file.read()

        last_index = 0

        for m in re.finditer(re_subroutine, data):

            subroutine_end = re.search(re_end_subroutine, data[m.end():]).start() + m.end()

            # assume comments are after the prototype
            comments = re.search(re_comments, data[m.end():subroutine_end])

            parameters = re.findall(re_parameters, data[m.end():subroutine_end])

            if comments and comments.start() <= 3:
                output_data += data[last_index:m.start()]
                if m.groups()[0] is not None:
                    subroutine_header = make_subroutine_header(comments.group(1), binds_to=m.group(1), parameters=parameters)
                else:
                    subroutine_header = make_subroutine_header(comments.group(1), parameters=parameters)

                output_data += subroutine_header
                output_data += data[m.start():m.end()+comments.start()]
                output_data += data[m.end()+comments.end():subroutine_end]

            else:
                output_data += data[last_index:m.start()]
                if m.groups()[0] is not None:
                    subroutine_header = make_subroutine_header(binds_to=m.group(1), parameters=parameters)
                else:
                    subroutine_header = make_subroutine_header(parameters=parameters)

                output_data += subroutine_header
                output_data += data[m.start():subroutine_end]

            last_index = subroutine_end


        output_data += data[last_index:]

    data = output_data
    output_data = ""
    last_index = 0

    re_function = re.compile(r"^[ \t]*function \S+\([\w ,&\n]*\)", flags=re.MULTILINE)

    re_end_function = re.compile(r"^[ \t]*end function", flags=re.MULTILINE)

    for m in re.finditer(re_function, data):

        function_end = re.search(re_end_function, data[m.end():]).start() + m.end()

        # assume comments are before the prototype
        comments = None
        for comments in re.finditer(re_comments,
                                    data[last_index:m.start() + 1]):
            pass

        parameters = re.findall(re_parameters, data[m.end():function_end])


        if comments and (m.start() - last_index - comments.end()) <= 3:
            output_data += data[last_index:last_index + comments.start()]
            function_header = make_function_header(comments.group(1), parameters=parameters)

        else:
            output_data += data[last_index:m.start()]
            function_header = make_function_header(parameters=parameters)

        last_index = function_end

        output_data += function_header
        output_data += data[m.start():function_end]

    output_data += data[last_index:]


    output_filename = filename + ".doxygen"

    with open(output_filename, 'w+') as output_file:
        output_file.write(output_data)



if __name__ == "__main__":
    filename = sys.argv[1]

    if filename[-2:] == ".H":
        process_header_file(filename)
    # elif filename[-4:] == ".cpp":
    #     process_cpp_file(filename)
    elif filename[-4:].lower() == ".f90":
        process_fortran_file(filename)
