import argparse

def main():
    # Create a parser
    parser = argparse.ArgumentParser(description="Process some integers.")

    # Add the "-fn" option
    parser.add_argument("-fn", "--filename", help="the name of the file")

    # Parse the command line arguments
    args = parser.parse_args()

    # Get the filename
    filename = args.filename

    # Now you can use the filename in your program
    print(f"The filename is {filename}")

if __name__ == "__main__":
    main()