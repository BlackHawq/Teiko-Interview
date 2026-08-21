.PHONY: setup pipeline dashboard

setup:
	python3 -m pip install -r requirements.txt

pipeline:
	python3 pipeline.py

dashboard:
	streamlit run dashboard.py
