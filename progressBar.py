from tqdm import tqdm


class ProgressBar:

    def __init__(self, total, description="Processing"):
        self.bar = tqdm(
            total=total,
            desc=description
        )

    def update(self, value):
        self.bar.update(value)

    def close(self):
        self.bar.close()