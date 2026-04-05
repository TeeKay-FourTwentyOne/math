"""
Generate 2000 candidate prompts for Experiment 0.
5 categories x 400 prompts each.
Output: prompts/candidates_2000.jsonl
"""
import json
import random
from pathlib import Path

SEED = 42
N_PER_CATEGORY = 400

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def L(s):
    """Parse comma-separated string into cleaned list."""
    return [x.strip() for x in s.strip().split(",") if x.strip()]


def ucfirst(s):
    """Uppercase first character, preserving rest."""
    return s[0].upper() + s[1:] if s else s


# ===========================================================================
# VALUE POOLS
# ===========================================================================

COUNTRIES = L("""
France, Japan, Brazil, Egypt, Australia, Germany, India, China, Russia, Mexico,
Canada, Italy, Spain, Portugal, Sweden, Norway, Finland, Denmark, Poland, Greece,
South Korea, Argentina, Peru, Colombia, Chile, Morocco, Tunisia, Kenya, Ghana,
Nigeria, Turkey, Iran, Pakistan, Thailand, Vietnam, Indonesia, Philippines,
Malaysia, Singapore, Nepal, Hungary, Czech Republic, Romania, Bulgaria, Croatia,
Ukraine, Kazakhstan, Mongolia, Cuba, Jamaica, Iceland, New Zealand, Sri Lanka,
Ethiopia, Tanzania, Bolivia, Ecuador, Venezuela, Lebanon, Jordan, Ireland
""")

BOOKS = L("""
Pride and Prejudice, 1984, To Kill a Mockingbird, The Great Gatsby,
One Hundred Years of Solitude, Moby Dick, War and Peace, Crime and Punishment,
The Catcher in the Rye, Brave New World, Animal Farm, Don Quixote,
Jane Eyre, Wuthering Heights, Great Expectations, The Count of Monte Cristo,
Les Miserables, Anna Karenina, Frankenstein, Dracula,
The Picture of Dorian Gray, Heart of Darkness, Huckleberry Finn,
Alice in Wonderland, Treasure Island, Robinson Crusoe, Gulliver's Travels,
A Tale of Two Cities, Oliver Twist, David Copperfield,
The Scarlet Letter, Little Women, The Three Musketeers, Ivanhoe,
The Metamorphosis, The Trial, The Stranger, The Plague,
Siddhartha, The Old Man and the Sea, Fahrenheit 451, Slaughterhouse-Five,
Lord of the Flies, The Grapes of Wrath, Of Mice and Men, Catch-22,
On the Road, Invisible Man, The Bell Jar, Beloved
""")

ELEMENTS = L("""
Hydrogen, Helium, Lithium, Beryllium, Carbon, Nitrogen, Oxygen, Fluorine,
Neon, Sodium, Magnesium, Aluminum, Silicon, Phosphorus, Sulfur, Chlorine,
Argon, Potassium, Calcium, Titanium, Chromium, Manganese, Iron, Cobalt,
Nickel, Copper, Zinc, Gallium, Germanium, Arsenic, Selenium, Bromine,
Krypton, Rubidium, Strontium, Silver, Tin, Iodine, Xenon, Barium,
Tungsten, Platinum, Gold, Mercury, Lead, Bismuth, Radon, Uranium,
Plutonium, Cesium
""")

EVENTS = L("""
French Revolution, American Civil War, World War I, World War II,
fall of the Berlin Wall, Apollo 11 Moon landing, Industrial Revolution,
Renaissance, Protestant Reformation, Russian Revolution,
signing of the Magna Carta, discovery of penicillin,
invention of the printing press, construction of the Panama Canal,
sinking of the Titanic, Great Fire of London, Black Death,
fall of the Roman Empire, first voyage of Columbus,
first powered flight by the Wright Brothers,
signing of the Declaration of Independence, Battle of Waterloo,
Cuban Missile Crisis, Chernobyl disaster, fall of Constantinople,
Boston Tea Party, Storming of the Bastille, Treaty of Versailles,
founding of the United Nations, launch of Sputnik,
assassination of Archduke Franz Ferdinand, Meiji Restoration,
partition of India, Korean War, Chinese Cultural Revolution,
Iranian Revolution, dissolution of the Soviet Union,
signing of the Treaty of Westphalia, Boxer Rebellion, Spanish Civil War,
Irish Potato Famine, abolition of apartheid, invention of the telephone,
discovery of the structure of DNA, first circumnavigation of the globe,
invention of the steam engine, signing of the Emancipation Proclamation,
eruption of Mount Vesuvius, discovery of radioactivity, Norman Conquest
""")

PEOPLE = L("""
Albert Einstein, Isaac Newton, Marie Curie, Charles Darwin, Nikola Tesla,
Leonardo da Vinci, William Shakespeare, Wolfgang Amadeus Mozart,
Ludwig van Beethoven, Mahatma Gandhi, Nelson Mandela,
Martin Luther King Jr., Abraham Lincoln, Cleopatra, Alexander the Great,
Julius Caesar, Napoleon Bonaparte, Galileo Galilei, Nicolaus Copernicus,
Johannes Kepler, Ada Lovelace, Florence Nightingale, Rosa Parks,
Frida Kahlo, Sigmund Freud, Karl Marx, Friedrich Nietzsche,
Aristotle, Plato, Socrates, Confucius, Genghis Khan, Marco Polo,
Thomas Edison, Alexander Graham Bell, Henry Ford, Alan Turing,
Grace Hopper, Pablo Picasso, Vincent van Gogh, Claude Monet, Rembrandt,
Johann Sebastian Bach, Pythagoras, Archimedes, Hippocrates, Sun Tzu,
Queen Victoria, Catherine the Great, Charlemagne
""")

ANIMALS = L("""
elephant, blue whale, cheetah, penguin, eagle, dolphin, octopus, gorilla,
tiger, polar bear, kangaroo, giraffe, wolf, fox, owl, shark, turtle,
salmon, butterfly, honeybee, ant, spider, chameleon, parrot, flamingo,
seahorse, jellyfish, starfish, koala, panda, sloth, armadillo, platypus,
hedgehog, rattlesnake, crocodile, hippopotamus, rhinoceros, zebra, ostrich,
hummingbird, peacock, lobster, squid, bat, moose, bison, lynx, otter, raven
""")

TOPICS_ABSTRACT = L("""
time, memory, silence, change, hope, fear, justice, truth, beauty,
loneliness, freedom, creativity, wisdom, patience, courage, empathy,
curiosity, ambition, nostalgia, resilience, gratitude, forgiveness,
purpose, identity, belonging, loss, growth, trust, wonder, simplicity,
chaos, harmony, doubt, faith, inspiration, solitude, adventure,
transformation, perseverance, imagination, consciousness, mortality,
infinity, symmetry, entropy, paradox, serendipity, melancholy,
euphoria, uncertainty
""")

TOPICS_CONCRETE = L("""
a sunrise over the desert, a thunderstorm at sea, an empty classroom at night,
a crowded subway car, a garden in autumn, a frozen lake at dawn,
a street market in summer, a library on a rainy day, a campfire under the stars,
a lighthouse in fog, an abandoned amusement park, a kitchen during a holiday meal,
a rooftop overlooking a city, a path through a dense forest, a train crossing a bridge,
a beach after a storm, a workshop filled with tools, a hospital waiting room,
a concert hall before the performance, an old bookshop, a harbor at sunset,
a mountain trail in spring, a cemetery in winter, a playground at dusk,
a bakery in the early morning, a river in flood, a field of sunflowers,
a cave with underground pools, a balcony during a rainstorm,
a greenhouse in January, a village after snowfall, a dock at low tide,
a windmill on a hill, a canyon at midday, a greenhouse full of orchids,
a cobblestone alley at midnight, a meadow full of fireflies,
a shipwreck on a reef, a clock tower at noon, a vineyard in October,
a waterfall hidden in jungle, a skating rink at night,
an attic full of old letters, a bridge over a misty river,
a volcano at twilight, a coral reef at depth, a sand dune in wind,
a cathedral interior, a barn during harvest, a telescope pointed at the sky,
a tide pool at low tide
""")

HYPOTHETICALS = L("""
you could travel to any point in history,
you woke up with the ability to read minds,
gravity suddenly became half as strong,
every person on Earth spoke the same language,
you could live underwater,
animals could talk to humans,
you never needed to sleep,
money ceased to exist tomorrow,
you could teleport anywhere instantly,
the internet disappeared permanently,
you had a perfect photographic memory,
humans could photosynthesize like plants,
you discovered a new habitable planet,
time moved at half its normal speed,
you could communicate with plants,
the Earth suddenly had two moons,
you could become invisible at will,
all borders between countries disappeared,
you could relive any single day of your life,
the seasons stopped changing forever,
you found a door to a parallel universe,
all humans could fly,
you could talk to your future self,
all music in the world was forgotten,
you were the last person on Earth,
you woke up a thousand years in the future,
rain started falling upward,
you could shrink to the size of an insect,
you discovered a message from aliens,
dreams became shared experiences,
color ceased to exist,
you could breathe in outer space,
shadows became solid objects,
you found a map to a hidden city,
the ocean turned to glass overnight,
books could write themselves,
trees could walk slowly,
silence became a visible substance,
you could pause time at will,
every lie became immediately obvious,
art came to life at night,
memories could be bottled and shared,
you could hear radio waves naturally,
snow fell warm instead of cold,
you found a staircase that went down forever,
you could taste colors,
the sun rose in the west,
you could split into two copies of yourself,
stones could speak, rocks remembered everything
""")

CREATIVE_TASKS = L("""
a holiday that celebrates failure, a sport played in zero gravity,
a language spoken only by children, a city built entirely underground,
a musical instrument made of ice, a school where animals teach humans,
a restaurant that serves food from the future, a library where books read to you,
a garden that grows in complete darkness, a vehicle powered by music,
a building that changes shape with the weather, a clock that measures emotions,
a map of an emotion, a recipe for courage, a machine that translates birdsong,
a sport played with shadows, a museum of forgotten smells,
a currency based on kindness, a calendar with thirteen months,
a game played only in dreams, a bridge connecting two memories,
a door that leads to yesterday, a window showing the view from any mountain,
a hat that lets you understand any animal, a book with a new ending each reading,
a song that changes the weather, a painting that moves when unobserved,
a key that opens any locked conversation, a mirror showing who you will be in ten years,
a blanket woven from starlight, an umbrella that shields you from sadness,
a compass that always points toward home, a telescope that sees into dreams,
a candle that burns brighter in darkness, a pen that only writes truth,
a chair that tells stories from its past, a tree growing different fruit each season,
a fountain granting one small wish per day, a mailbox receiving letters from the past,
a staircase leading to a different era, a bell that rings for approaching danger,
a net that catches forgotten words, a ship that sails through clouds,
a train that travels through paintings, a bookmark that pulls you into the story,
a pillow that records your dreams, a lantern powered by laughter,
a roof that opens for meteor showers, a door knocker announcing visitors by intent,
a well that echoes future conversations
""")

# Pools for continuation prompts
ADJECTIVES = L("""
old, quiet, abandoned, bustling, mysterious, forgotten, ancient, remote,
crumbling, sunlit, fog-shrouded, windswept, isolated, crowded, elegant,
modest, sprawling, hidden, weathered, charming, derelict, ornate, narrow,
ivy-covered, dimly lit, whitewashed, moss-covered, painted, locked, towering
""")

NOUN_PLACES = L("""
lighthouse, library, train station, cafe, hospital, warehouse, cathedral,
farmhouse, observatory, courthouse, theater, museum, market, harbor,
laboratory, university, monastery, factory, bridge, tavern,
chapel, bookshop, clocktower, inn, manor, cottage, gallery, forge,
mill, schoolhouse
""")

CHARACTERS = L("""
a young woman, an elderly professor, a traveling merchant, a retired soldier,
a curious child, a nervous clerk, a foreign diplomat, a local fisherman,
a lone traveler, a street musician, a night watchman, a retired doctor,
a school teacher, a wandering monk, a newspaper reporter, a ship captain,
a clockmaker, a librarian, a blacksmith, a painter,
a detective, a botanist, a cartographer, a photographer, a translator,
a pilot, a baker, an architect, a jeweler, a violinist
""")

TIMES_OF_DAY = L("""
early morning, mid-morning, noon, early afternoon, late afternoon,
dusk, early evening, late evening, midnight, the small hours of the morning,
just before dawn, twilight, the golden hour, midday, nightfall
""")

SETTINGS = L("""
a small coastal town, the university campus, the financial district,
a mountain village, the industrial quarter, a suburban neighborhood,
an island community, the old quarter of the city, a farming valley,
a border town, a lakeside resort, a mining settlement,
the artists' district, a fishing port, a desert outpost,
a riverside market town, a hill station, a frontier settlement,
a monastery perched on sea cliffs, a city built on canals
""")

OBJECTS = L("""
a leather-bound journal, a brass compass, a skeleton key, an old photograph,
a sealed envelope, a hand-drawn map, a wooden music box, a silver locket,
a glass vial, a carved figurine, a faded letter, a pocket watch,
a small painting, a wax-sealed document, a crystal paperweight,
a tarnished mirror, a set of ivory dice, a silk handkerchief,
a rusted lockbox, a rolled parchment
""")

PROFESSIONS = L("""
detective, archaeologist, marine biologist, meteorologist, historian,
journalist, cryptographer, surgeon, geologist, linguist,
astronomer, folklorist, epidemiologist, civil engineer, archivist
""")

PHENOMENA = L("""
strange lights in the sky, unusual sounds after midnight,
missing pets in the neighborhood, unexplained marks on walls,
recurring power outages, anonymous letters appearing in mailboxes,
sudden temperature drops in certain rooms, shadows moving independently,
clocks stopping at the same time, birds flying in unnatural patterns,
an unfamiliar scent drifting through town, objects rearranging overnight,
footprints appearing in fresh snow, radio signals from unknown sources,
plants growing at unnatural speed, water levels changing without rain,
animals behaving strangely, dust collecting in geometric patterns,
echoes in rooms that should not have them, windows fogging without cause
""")

# Pools for conversational prompts
CASUAL_TOPICS = L("""
pineapple on pizza, working from home, social media, artificial intelligence,
climate change, space exploration, electric cars, online education,
vegetarianism, universal basic income, self-driving cars, cryptocurrency,
fast fashion, mental health awareness, remote work, genetic engineering,
nuclear energy, video games as art, standardized testing,
the four-day work week, pet ownership in cities, public transportation,
streaming services, handwriting versus typing, breakfast for dinner,
small talk with strangers, tipping culture, daylight saving time,
reality television, audiobooks versus physical books, morning routines,
minimalism, cooking at home, professional sports salaries,
music festivals, urban gardening, tiny houses, board games as a hobby,
learning a new language as an adult, solo travel,
vintage fashion, specialty coffee, digital detox, volunteering,
thrift shopping, meditation, journaling, outdoor fitness,
podcasts, food trucks
""")

HELP_TASKS = L("""
write a thank-you note to my mentor, plan a birthday party on a budget,
organize my home office, choose a laptop for college,
prepare for a job interview, start a vegetable garden,
learn basic car maintenance, write a cover letter for a new job,
create a weekly meal plan, set up a home network, plan a road trip,
choose a good book for my book club, write a best man speech,
pack efficiently for a two-week trip, start a morning exercise routine,
organize a community cleanup event, learn to cook Italian food,
set up a personal budget spreadsheet, write a letter of recommendation,
plan a surprise anniversary dinner, choose the right running shoes,
prepare for my first marathon, learn basic photography skills,
write a short story for a contest, plan a small home renovation,
start investing with a small budget, find a volunteer opportunity,
learn to play acoustic guitar, set up a freshwater aquarium,
create a study schedule for exams, plan a weekend camping trip,
choose the right breed of dog to adopt, write a formal complaint letter,
learn calligraphy basics, plan a dinner party for eight people,
fix a leaky kitchen faucet, build a simple bookshelf,
choose a paint color for my living room, start a personal blog,
learn basic sewing skills, plan a garden layout for spring,
choose a good day-hike trail, write meaningful wedding vows,
prepare a ten-minute presentation, organize a neighborhood garage sale,
learn origami for beginners, create a vacation itinerary,
start a home composting system, set up a cozy reading nook,
learn to play a new card game
""")

RECOMMEND_ITEMS = L("""
book to read this winter, podcast about science, documentary about nature,
recipe for a beginner cook, workout routine for busy people,
app for learning a new language, museum to visit in Europe,
board game for a group of adults, houseplant for a low-light apartment,
online course about world history, movie for a rainy weekend afternoon,
trail for a beginner day hike, journal for daily reflective writing,
tool for a home improvement project, thoughtful gift for a book lover,
dish to try at a Thai restaurant, hobby to pick up for stress relief,
travel destination in Southeast Asia, first programming language to learn,
song to listen to while studying, coffee brewing method for home,
way to meet people when new to a city, useful skill to learn in a weekend,
healthy and quick breakfast recipe, true crime podcast series,
strategy for getting better sleep, exercise that helps with back pain,
accessible philosophy book, approach to practicing mindfulness,
game for a long car ride, art supply set for a complete beginner,
resource for learning chess strategy, city to visit in South America,
technique for improving vocabulary, high-energy snack for hiking,
show to binge-watch this weekend, effective charity to donate to,
instrument to learn as an adult, type of herbal tea to try,
way to decorate a very small apartment, cultural festival worth attending,
relaxing craft project for weekends, digital tool for taking better notes,
method for reading faster, low-impact sport to start as a beginner,
way to reduce daily screen time, meal prep recipe for the work week,
fun approach to learning math, magazine worth subscribing to,
documentary about space exploration
""")


# ===========================================================================
# GENERATORS
# ===========================================================================

def gen_factual():
    r = random.Random(SEED)
    p = []
    # 8 template groups, ~50 each
    for c in r.sample(COUNTRIES, 50):
        p.append(f"What is the capital of {c}?")
    for b in r.sample(BOOKS, 50):
        p.append(f"Who wrote the novel '{b}'?")
    for e in r.sample(ELEMENTS, 50):
        p.append(f"What is the chemical symbol for {e}?")
    for ev in r.sample(EVENTS, 50):
        p.append(f"In what year did the {ev} occur?")
    for pe in r.sample(PEOPLE, 50):
        p.append(f"What is {pe} known for?")
    for c in r.sample(COUNTRIES, 50):
        p.append(f"What language is primarily spoken in {c}?")
    for c in r.sample(COUNTRIES, 50):
        p.append(f"What continent is {c} located on?")
    for a in r.sample(ANIMALS, 50):
        p.append(f"What is the typical lifespan of a {a}?")
    p = list(dict.fromkeys(p))
    r.shuffle(p)
    return [{"text": t, "category": "factual"} for t in p[:N_PER_CATEGORY]]


def gen_open_ended():
    r = random.Random(SEED + 1)
    p = []
    for t in r.sample(TOPICS_ABSTRACT, 50):
        p.append(f"Write a paragraph about {t}.")
    for t in r.sample(TOPICS_CONCRETE, 50):
        p.append(f"Describe {t} in vivid detail.")
    for h in r.sample(HYPOTHETICALS, 50):
        p.append(f"Imagine {h}. What would your first day be like?")
    for c in r.sample(CREATIVE_TASKS, 50):
        p.append(f"Invent and describe {c}.")
    for t in r.sample(TOPICS_ABSTRACT, 50):
        p.append(f"How would you explain {t} to a five-year-old?")
    for a in r.sample(ANIMALS, 50):
        p.append(f"Write a short story from the perspective of a {a}.")
    # Compare pairs (filter self-comparisons)
    a_list = r.sample(TOPICS_ABSTRACT, len(TOPICS_ABSTRACT))
    b_list = r.sample(TOPICS_ABSTRACT, len(TOPICS_ABSTRACT))
    for a, b in zip(a_list, b_list):
        if a != b:
            p.append(f"Compare and contrast {a} and {b}.")
    for h in r.sample(HYPOTHETICALS, 50):
        p.append(f"What would happen if {h}? Describe the consequences.")
    p = list(dict.fromkeys(p))
    r.shuffle(p)
    return [{"text": t, "category": "open_ended"} for t in p[:N_PER_CATEGORY]]


def gen_reasoning():
    r = random.Random(SEED + 2)
    p = []

    # Addition (55)
    for _ in range(55):
        a, b = r.randint(10, 9999), r.randint(10, 9999)
        p.append(f"What is {a} + {b}?")

    # Subtraction (55)
    for _ in range(55):
        a = r.randint(100, 9999)
        b = r.randint(10, a)
        p.append(f"What is {a} - {b}?")

    # Multiplication (50)
    for _ in range(50):
        a, b = r.randint(2, 99), r.randint(2, 99)
        p.append(f"What is {a} times {b}?")

    # Division (25)
    for _ in range(25):
        b = r.randint(2, 50)
        a = b * r.randint(2, 200)
        p.append(f"What is {a} divided by {b}?")

    # Percentage (25)
    for _ in range(25):
        pct = r.choice([5, 10, 15, 20, 25, 30, 40, 50, 60, 75, 80, 90])
        n = r.randint(50, 5000)
        p.append(f"What is {pct}% of {n}?")

    # Distance-speed-time (50)
    for _ in range(17):
        speed = r.randint(20, 120)
        hours = r.randint(1, 12)
        p.append(f"A car travels at {speed} miles per hour for {hours} hours. How far does it go?")
    for _ in range(17):
        dist = r.randint(50, 500)
        speed = r.randint(20, 100)
        p.append(f"How many hours does it take to travel {dist} miles at {speed} miles per hour?")
    for _ in range(16):
        dist = r.randint(100, 1000)
        hours = r.randint(2, 10)
        p.append(f"A train covers {dist} miles in {hours} hours. What is its average speed?")

    # Money (50)
    for _ in range(17):
        items = r.randint(2, 20)
        price = round(r.uniform(0.50, 50.0), 2)
        p.append(f"You buy {items} items at ${price:.2f} each. What is the total cost?")
    for _ in range(17):
        budget = r.randint(50, 500)
        price = round(r.uniform(5.0, 50.0), 2)
        p.append(f"You have ${budget}. Each item costs ${price:.2f}. How many can you buy?")
    for _ in range(16):
        original = r.randint(20, 200)
        discount = r.choice([10, 15, 20, 25, 30, 40, 50])
        p.append(f"A shirt originally costs ${original}. It is {discount}% off. What is the sale price?")

    # Sequences (50)
    for _ in range(25):
        start = r.randint(1, 50)
        step = r.randint(1, 15)
        seq = [start + i * step for i in range(5)]
        p.append(f"What comes next in the sequence: {', '.join(map(str, seq))}, ?")
    for _ in range(25):
        start = r.randint(1, 10)
        ratio = r.randint(2, 5)
        seq = [start * ratio ** i for i in range(5)]
        p.append(f"What comes next in the sequence: {', '.join(map(str, seq))}, ?")

    # Multi-step (50)
    for _ in range(25):
        people_orig = r.randint(4, 12)
        people_new = r.randint(people_orig + 1, people_orig * 3)
        cups = r.randint(2, 8)
        p.append(
            f"A recipe calls for {cups} cups of flour to serve {people_orig} people. "
            f"How much flour do you need to serve {people_new} people?"
        )
    for _ in range(25):
        n_workers = r.randint(3, 10)
        days = r.randint(5, 30)
        new_workers = r.randint(n_workers + 1, n_workers * 3)
        p.append(
            f"If {n_workers} workers can finish a job in {days} days, "
            f"how many days would {new_workers} workers take?"
        )

    p = list(dict.fromkeys(p))
    r.shuffle(p)
    return [{"text": t, "category": "reasoning"} for t in p[:N_PER_CATEGORY]]


def gen_continuation():
    r = random.Random(SEED + 3)
    p = []

    # Template 1: Discovery (50)
    for _ in range(50):
        adj = r.choice(ADJECTIVES)
        place = r.choice(NOUN_PLACES)
        char = r.choice(CHARACTERS)
        p.append(
            f"The {adj} {place} had been locked for as long as anyone could remember. "
            f"When {char} finally found a way inside,"
        )

    # Template 2: Time and setting (50)
    for _ in range(50):
        tod = r.choice(TIMES_OF_DAY)
        setting = r.choice(SETTINGS)
        char = r.choice(CHARACTERS)
        p.append(
            f"It was {tod} in {setting}. "
            f"{ucfirst(char)} was walking alone when something caught their eye:"
        )

    # Template 3: Object with history (50)
    for _ in range(50):
        obj = r.choice(OBJECTS)
        adj = r.choice(ADJECTIVES)
        place = r.choice(NOUN_PLACES)
        char = r.choice(CHARACTERS)
        p.append(
            f"{ucfirst(obj)} had been sitting in the {adj} {place} for years, untouched. "
            f"When {char} finally picked it up,"
        )

    # Template 4: Mystery (50)
    for _ in range(50):
        char = r.choice(CHARACTERS)
        adj = r.choice(ADJECTIVES)
        place = r.choice(NOUN_PLACES)
        p.append(
            f"Nobody could explain why {char} kept returning to the {adj} {place}. "
            f"Some said it was because"
        )

    # Template 5: Aftermath (50)
    for _ in range(50):
        place = r.choice(NOUN_PLACES)
        setting = r.choice(SETTINGS)
        char = r.choice(CHARACTERS)
        p.append(
            f"Three days after the incident at the {place}, {setting} was still on edge. "
            f"{ucfirst(char)} had been asking questions that nobody wanted to answer."
        )

    # Template 6: Unexpected find (50)
    for _ in range(50):
        char = r.choice(CHARACTERS)
        adj = r.choice(ADJECTIVES)
        place = r.choice(NOUN_PLACES)
        obj = r.choice(OBJECTS)
        p.append(
            f"The last thing {char} expected to find in the {adj} {place} was {obj}. "
            f"And yet there it was, sitting in plain sight."
        )

    # Template 7: Investigation (50)
    for _ in range(50):
        phenom = r.choice(PHENOMENA)
        prof = r.choice(PROFESSIONS)
        setting = r.choice(SETTINGS)
        p.append(
            f"Reports of {phenom} had been increasing for weeks in {setting}. "
            f"A {prof} was dispatched to investigate. Upon arrival,"
        )

    # Template 8: Dialogue (50)
    for _ in range(50):
        char = r.choice(CHARACTERS)
        adj = r.choice(ADJECTIVES)
        place = r.choice(NOUN_PLACES)
        p.append(
            f"\"I wouldn't go near the {adj} {place} if I were you,\" "
            f"said {char}, lowering their voice. \""
        )

    p = list(dict.fromkeys(p))
    r.shuffle(p)
    return [{"text": t, "category": "continuation"} for t in p[:N_PER_CATEGORY]]


def gen_conversational():
    r = random.Random(SEED + 4)
    p = []

    # Opinion (60)
    for t in r.sample(CASUAL_TOPICS, 50):
        p.append(f"What do you think about {t}?")
    for t in r.sample(CASUAL_TOPICS, 10):
        p.append(f"I'm curious -- what's your take on {t}?")

    # Help (60)
    for t in r.sample(HELP_TASKS, 50):
        p.append(f"Can you help me {t}?")
    for t in r.sample(HELP_TASKS, 10):
        p.append(f"I could really use some help. I need to {t}.")

    # Recommendations (50)
    for t in r.sample(RECOMMEND_ITEMS, 50):
        p.append(f"Can you recommend a good {t}?")

    # Explanations (50)
    for t in r.sample(CASUAL_TOPICS, 50):
        p.append(f"Can you explain {t} in simple terms?")

    # Greetings (50)
    greetings = ["Hey", "Hi there", "Hello", "Good morning", "Good evening",
                 "Howdy", "Hey there", "Hi", "Greetings", "Good afternoon"]
    follow_ups = [
        "How's it going?", "What's new?", "I have a question for you.",
        "Can I ask you something?", "I need your help with something.",
    ]
    for g in greetings:
        for f in follow_ups:
            p.append(f"{g}! {f}")

    # Hypothetical conversations (50)
    for h in r.sample(HYPOTHETICALS, 50):
        p.append(f"Just a fun question -- if {h}, what would you do?")

    # Advice (40)
    advice_situations = [
        "changing careers in my thirties",
        "dealing with a difficult coworker",
        "managing my time better as a student",
        "staying motivated while working from home",
        "handling criticism without taking it personally",
        "saving money on a tight budget",
        "making friends in a new city",
        "overcoming public speaking anxiety",
        "balancing work and personal life",
        "deciding whether to go back to school",
        "learning to say no without feeling guilty",
        "dealing with burnout at work",
        "building better daily habits",
        "improving my writing skills",
        "having a difficult conversation with a family member",
        "reducing stress during exam season",
        "choosing between two good job offers",
        "getting better at cooking for one person",
        "starting an exercise routine and sticking to it",
        "dealing with imposter syndrome",
        "improving my sleep schedule",
        "networking when you are introverted",
        "managing finances as a freelancer",
        "learning to be more patient",
        "setting boundaries in relationships",
        "transitioning to a plant-based diet",
        "handling a long-distance relationship",
        "picking up a creative hobby",
        "dealing with information overload",
        "making big life decisions with uncertainty",
        "getting organized when you tend toward chaos",
        "navigating a career change into tech",
        "being more present and less distracted",
        "building self-confidence gradually",
        "coping with seasonal mood changes",
        "learning to enjoy being alone",
        "managing anxiety about the future",
        "developing a growth mindset",
        "finding meaning in routine work",
        "becoming a better listener",
    ]
    for s in advice_situations:
        p.append(f"I need some advice about {s}. What would you suggest?")

    # Reflective (40)
    for t in r.sample(TOPICS_ABSTRACT, 40):
        p.append(f"What do you think is the most important thing about {t}?")

    p = list(dict.fromkeys(p))
    r.shuffle(p)
    return [{"text": t, "category": "conversational"} for t in p[:N_PER_CATEGORY]]


# ===========================================================================
# MAIN
# ===========================================================================

def main():
    out_dir = Path(__file__).parent / "prompts"
    out_dir.mkdir(parents=True, exist_ok=True)

    all_prompts = []
    generators = [
        ("factual",       gen_factual),
        ("open_ended",    gen_open_ended),
        ("reasoning",     gen_reasoning),
        ("continuation",  gen_continuation),
        ("conversational", gen_conversational),
    ]
    for name, gen_fn in generators:
        batch = gen_fn()
        print(f"  {name}: {len(batch)} prompts")
        all_prompts.extend(batch)

    # Assign global indices
    for i, p in enumerate(all_prompts):
        p["idx"] = i

    # Dedup sanity check
    texts = [p["text"] for p in all_prompts]
    unique = len(set(texts))
    if unique < len(texts):
        print(f"WARNING: {len(texts) - unique} duplicate prompts (will keep first occurrence)")
        seen = set()
        deduped = []
        for p in all_prompts:
            if p["text"] not in seen:
                seen.add(p["text"])
                deduped.append(p)
        all_prompts = deduped
        for i, p in enumerate(all_prompts):
            p["idx"] = i

    # Save
    out_path = out_dir / "candidates_2000.jsonl"
    with open(out_path, "w", encoding="utf-8") as f:
        for p in all_prompts:
            f.write(json.dumps(p, ensure_ascii=False) + "\n")

    from collections import Counter
    cats = Counter(p["category"] for p in all_prompts)
    print(f"\nGenerated {len(all_prompts)} prompts total:")
    for cat, count in sorted(cats.items()):
        print(f"  {cat}: {count}")
    print(f"\nSaved to {out_path}")

    # Show a few samples per category
    print("\nSamples:")
    r = random.Random(SEED)
    for cat in sorted(cats.keys()):
        cat_prompts = [p for p in all_prompts if p["category"] == cat]
        samples = r.sample(cat_prompts, min(2, len(cat_prompts)))
        for s in samples:
            text = s["text"][:90] + ("..." if len(s["text"]) > 90 else "")
            print(f"  [{cat}] {text}")


if __name__ == "__main__":
    main()
