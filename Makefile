.PHONY: setup pipeline dashboard

setup:
	pip install -r requirements.txt

pipeline:
	python load_data.py
	python analysis/part2_frequencies.py
	python analysis/part3_stats.py
	python analysis/part4_subset.py

dashboard:
	streamlit run dashboard/app.py
