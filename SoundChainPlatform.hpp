#pragma once

#include "SoundChain.hpp"

struct SoundChainPlatformSettings {
	int SampleRate;
	int Channels;
};

// TODO: add buffering (double/triple or circular)

class SoundChainPlatform {
protected:
	SoundChainPlatformSettings _settings;
	unsigned int _samplesElapsed;

	virtual void Setup() {};
	virtual void Start() {};
	virtual void End() {};

public:
	SoundChainPlatform() {}
	~SoundChainPlatform() {}

	void Initialize(SoundChainPlatformSettings soundChainPlatformSettings) {
		if (_initialized) return;

		_settings = soundChainPlatformSettings;

		Setup();

		if (_previous) {
			_previous->Initialize(GetSoundChainSettings());
		}

		Start();

		_samplesElapsed = 0.0f;
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

	SoundChainSettings GetSoundChainSettings() {
		SoundChainSettings settings = {_settings.SampleRate, _settings.Channels};
		return settings;
	}

	unsigned int samplesElapsed() {return _samplesElapsed;}
	double time() {return (double)_samplesElapsed / (double)_settings.SampleRate;}

private:
	SoundChainBase* _previous = nullptr;
	bool _initialized = false;
};