from tqdm.auto import tqdm

class ProcessTracker:
    """Responsible for tracking and displaying progess """
    def __init__(self, total_items, description="Processing"):
        self.bar = tqdm(total=total_items, desc=description)
    
    def update(self, increment=1):
        self.bar.update(increment)

    def close(self):
        self.bar.close()


