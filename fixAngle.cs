using System;
using UnityEngine;



/// <summary>
/// A representtaation of an angle, always under 360 degrees and with a easy toRad function.
/// </summary>
public class FixAngle
{
	private float _angle;
	
	public FixAngle(float angle)
	{
		_angle = angle % 360;
	}
	
	public float angle()
	{
		return _angle;
	}
	
	public float toRad(){
		return _angle * Mathf.Deg2Rad;
	}
}

