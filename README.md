# protostar-variability


NOTES: 
epochs table over time: 
- 29-06-23 epochs table:
    - this splits epochs both by date (every 15 days) and makes a new epoch if n_images>100
    - not ideal bc some epochs only have a few images in them after this arbitrary cut
    - So: I want to go back to the old version, and handle cases individually in the big loop
 
- previous fixes:
    - L483 coordinates weren't centered
    - 
