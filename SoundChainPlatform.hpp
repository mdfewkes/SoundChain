#pragma once

#include "SoundChain.hpp"

class SoundChainPlatform {
protected:
	SoundChainSettings _settings;
	unsigned int _samplesElapsed;

	virtual void Setup() {};
	virtual void Start() {};
	virtual void End() {};

public:
	SoundChainPlatform() {}
	~SoundChainPlatform() {}

	void Initialize(SoundChainSettings &soundChainSettings) {
		if (_initialized) return;

		_settings = soundChainSettings;

		Setup();

		if (_previous) {
			_previous->Initialize(_settings);
		}

		Start();

		_initialized = true;
	}

	bool IsInitialized() {return _initialized;}

	void Terminate() {
		if (!_initialized) return;

		End();

		_initialized = false;
	}

	SoundChainBase* GetPrevious() {return _previous;}
	void SetPrevious(SoundChainBase* previous) {_previous = previous;}

	unsigned int sampleElapsed() {return _samplesElapsed;}
	double time() {return (double)_samplesElapsed / (double)_settings.SampleRate;}

private:
	SoundChainBase* _previous = nullptr;
	bool _initialized = false;
};