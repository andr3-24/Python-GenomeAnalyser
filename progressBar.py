from tqdm import tqdm


class ProgressBar:

    def __init__(self, total, description="Loading.."):

        """
        Initialize the progress bar object.

        Parameters:
        total (int): The total number of steps/items that need to be processed.
        description (str): The text displayed before the progress bar.
        """
        
        
        '''
        self is a reference to the current object instance.
        It allows the class methods to access and modify the attributes
        that belong to this specific object.
        
        In this case, self.bar stores the tqdm progress bar created
        during object initialization, so other methods like update()
        and close() can use the same progress bar instance.
        '''

        # Create a tqdm progress bar instance and store it as an attribute
        # of the ProgressBar object.
        #
        # self.bar allows other methods of this class (update, close)
        # to access and modify the same progress bar.
        self.bar = tqdm(
            total=total,
            desc=description,
            unit="step",
            colour="GREEN"
        )


    def update(self, amount=1):
        """
        Update the progress bar.

        Parameters:
        amount (int): Number of completed steps to add to the progress bar.
        """

        # Increase the progress bar progress by the given amount.
        # Default value is 1, meaning one step is completed.
        self.bar.update(amount)


    def close(self):
        """
        Close the progress bar.
        """

        # Properly close the tqdm progress bar.
        # This releases resources and leaves the terminal in a clean state.
        self.bar.close()